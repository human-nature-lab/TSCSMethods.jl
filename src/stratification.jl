# stratified.jl

#=
all of the below essentially follow from splitting matches

just execute on each group of a grouped dataframe,
and append output
=#

"""
    stratify!(stratfunc::Function, args...; kwargs...)

Apply a stratification function, its arguments, to apply the stratification, calculate the stratified grandbalances, the treated observations in each group, and return the updated model with plot labels.
"""
function stratify(
  stratfunc::Function, model::CIC, args...; kwargs...
)

  strata, stratlabels, stratifier = stratfunc(model, args...; kwargs...)

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, matches, balances, meanbalances, iterations, estimator = model;

  stratmodel = CICStratified(
    title = title,
    id = id,
    t = t,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    timevary = timevary,
    stratifier = stratifier,
    strata = strata,
    reference = reference,
    F = F, L = L,
    observations = observations,
    ids = ids,
    matches = matches,
    balances = balances,
    meanbalances = meanbalances,
    grandbalances = GrandDictStrat(),
    iterations = iterations,
    results = DataFrame(),
    treatednum = Dict{Int64, Int64}(), ##
    estimator = estimator,
    labels = stratlabels
  );
  
  if nrow(meanbalances) == 0
    error("calculate meanbalances first")
  end

  grandbalance!(stratmodel);
  treatednums!(stratmodel);

  return stratmodel
end

"""
assign treatment numbers in each category
"""
function treatednums!(model)
  @unpack strata, treatednum = model;
  for s in unique(strata)
    treatednum[s] = sum(strata .== s)
  end
  return model
end

"""
    customstrat!(
      cc,
      stratdict::Union{Dict{Tuple{Int64, Int64}, Int64}, Dict{Int64, Int64}}
    )

Stratify based on the values of some input dictionary, specifying strata for each (t, id) or each (id) for a stratification that varies only by unit.
"""
function customstrat(
  model, stratifier,
  stratdict::Dict{Tuple{Int64, Int64}, Int64}
)

  @unpack observations = model;

  for (i, ob) in enumerate(observations)
    strata[i] = get(stratdict, ob, 0)
  end

  return strata, Dict{Int, String}(), stratifier # labels should come from elsewhere
end

function assignq(val, X, lenX)
  q = 0
  for i in 1:(lenX-1)
    if val >= X[i]
      q += 1
    end
  end
  return q
end

"""
Generic function to stratify on a covariate present in the dataframe.

Assumes that matching, balancing, and meanbalancing has ocurred.
"""
function variablestrat(
  model::CIC, stratifier, dat;
  qtes = [0, 0.25, 0.5, 0.75, 1.0], zerosep = false
)

  @unpack title, id, t, outcome, treatment, timevary, observations = model;

  strata = Vector{Int}(undef, length(observations));

  c1 = dat[:, treatment] .== 1;

  missingpresent = false;

  # if separate zero
  zs = 0;

  if !timevary[stratifier]
    udf = unique(@view(dat[c1, :]), [id, stratifier], view = true);
        
    udict = Dict{Tuple{Int64, Int64}, eltype(udf[!, stratifier])}();
    @eachrow udf begin
      udict[($(t), $(id))] = $stratifier
    end
    
    xvec = udf[!, stratifier];
    missingpresent = any(ismissing.(udf[!, stratifier])) # update value
    # if there are missing values in the data, just skip them
    if missingpresent
      xvec = disallowmissing(xvec[.!ismissing.(xvec)]);
    end
    X = sort(quantile(xvec, qtes));
    Xlen = length(X);

    for (i, ob) in enumerate(observations)
      obval = get(udict, ob, 0)
      strata[i] = !ismissing(obval) ? assignq(obval, X, Xlen) : Xlen
    end
    
  elseif timevary[stratifier] # do at time of treatment
    
    X = sort(quantile(@views(dat[c1, var])));
    Xlen = length(X);

    udf = unique(@view(dat[c1, :]), [id, stratifier], view = true);

    # treatment events
    udict = Dict{Tuple{Int64, Int64}, eltype(udf[!, stratifier])}();
    @eachrow udf begin
      udict[($(t), $(id))] = $stratifier
    end
    
    if zerosep
      c2 = dat[:, var] .> 0;
      X = sort(quantile(@views(dat[c2 .& c1, var])));
      zs = 1;
    end

    for (i, ob) in enumerate(observations)
      obval = get(udict, ob, 0)
      strata[i] = !ismissing(obval) ? assignq(obval, X, Xlen) : Xlen
    end
    
  end

  stratlabels = label_variablestrat(
    string.(round.(X, digits = 2));
    missingpresent = missingpresent,
    zerosep = zerosep
  )

  return strata, stratlabels, stratifier
end

function label_variablestrat(
  quantnames; missingpresent = false, zerosep = false
)

  zs = !zerosep ? 0 : 1

  labels = Dict{Int, String}();
  sizehint!(labels, length(quantnames) - 1);
  rnge = (1:length(quantnames) - 1) .+ zs
  for i in rnge
    labels[i] = quantnames[i - zs] * " to " * quantnames[i + 1 - zs]
  end

  if zerosep
    labels[1] = "0.0"
  end

  if missingpresent
    labels[length(quantnames) + zs] = "Missing Values"
  end
  return labels
end

function combostratname(vars)
  varn = ""
  for (i, var) in enumerate(vars)
    if i < length(vars)
      varn = varn * string(var) * " x "
    else
      varn = varn * string(var)
    end
  end
  
  return Symbol(varn)
end

"""
    combostrat!(cc, dat, vars::Vector{Symbol}; varslabs = nothing)

Stratify based on the combinations of one or more variables. Strata are formed directly from the variable values.
"""
function combostrat(cc, dat, vars::Vector{Symbol}; varslabs = nothing)

  ### example vars
  # dat[!, :hightrump] = dat[!, vn.ts16] .>= 0.50;
  # dat[!, :highinc] = dat[!, vn.mil] .>= tscsmethods.mean(dat[!, vn.mil]);

  # vars = [:hightrump, :highinc]
  ###
  c1 = dat[:, cc.treatment] .== 1;
  udf = unique(@view(dat[c1, :]), [cc.t, cc.id, vars...], view = true);
  
  combos = Iterators.product([unique(udf[!, var]) for var in vars]...);
  
  stratmap = Dict(
    reshape(collect(combos), length(combos)) .=> 1:length(combos)
  );
  
  udict = Dict{Tuple{Int, Int}, Int}();
  for r in eachrow(udf)
    udict[(r[cc.t], r[cc.id])] = stratmap[(r[vars]...)]
  end
  
  for (i, ob) in enumerate(observations)
    obval = get(udict, ob, 0)
    strata[i] = !ismissing(obval) ? assignq(obval, X, Xlen) : Xlen
  end
  
  stratathatexist = unique(cc.meanbalances.stratum);

  stratlabels = Dict{Int, String}();
  for (k, v) in stratmap
    if stratmap[k] ∈ stratathatexist
      stratlabels[v] = combostratlab(vars, k, varslabs)
    end
  end
  
  return strata, stratlabels. combostratname(vars)
end

function combostratlab(vars, k, varslabs)
  init = ""

  for i in eachindex(vars)
    varlab = get(varslabs, vars[i], missing)
    vi = string(vars[i])
    if !ismissing(varlab)
      ki = string(get(varlab, k[i], missing))
    else
      ki = string(k[i])
    end
    
    if i < length(vars)
      init = init * vi * " = " * ki * ", "
    else
      init = init * vi * " = " * ki
    end
  end
  return init
end

function combostratlab(vars, k)
  init = ""
  for i in eachindex(vars)
    if i < length(vars)
      init = init * string(vars[i]) * " = " * string(k[i]) * ", "
    else
      init = init * string(vars[i]) * " = " * string(k[i])
    end
  end
  return init
end
