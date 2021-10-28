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
function stratify!(stratfunc::Function, args...; kwargs...)

  cc, labels = stratfunc(args...; kwargs...)
  
  grandbalance!(cc)
  treatednums!(cc)

  return cc, labels
end

"""
assign treatment numbers in each category
"""
function treatednums!(cc::cicmodel)
  cc.treatednum = Dict{Int64, Int64}();
  for s in unique(cc.meanbalances.stratum)
    
    cc.treatednum[s] = nrow(
      unique(
        @view cc.meanbalances[!, [:treattime, :treatunit]][cc.meanbalances.stratum .== s, :]
      )
    )
  end
  return cc
end

"""
assign treatment numbers left in each category for caliper model
"""
function treatedlefts!(cc::AbstractCICModel)
  cc.treatedleft = Dict{Int64, Int64}();
  for s in unique(cc.meanbalances.stratum)
    
    cc.treatedleft[s] = nrow(
      unique(@view cc.meanbalances[!, [:treattime, :treatunit]][cc.meanbalances.stratum .== s, :])
    )
  end
  return cc
end

# generic stratification functions

"""
    customstrat!(
      cc,
      stratdict::Union{Dict{Tuple{Int64, Int64}, Int64}, Dict{Int64, Int64}}
    )

Stratify based on the values of some input dictionary, specifying strata for each (t, id) or each (id) for a stratification that varies only by unit.
"""
function customstrat!(
  cc, stratifier,
  stratdict::Union{Dict{Tuple{Int64, Int64}, Int64}, Dict{Int64, Int64}}
)

  if typeof(stratdict) == Dict{Tuple{Int64, Int64}, Int64}
    cc.matches[!, :stratum] = Vector{Int}(undef, nrow(cc.matches));
    @eachrow! cc.matches begin
      :stratum = stratdict[(:treattime, :treatunit)]
    end

    cc.meanbalances[!, :stratum] = Vector{Int}(undef, nrow(cc.meanbalances));
    @eachrow! cc.meanbalances begin
    :stratum = stratdict[(:treattime, :treatunit)]
    end
  else
    @eachrow! cc.matches begin
      # includes zerosep case by loop design (q = 0)
      :stratum = stratdict[:treatunit]
    end

    # meanbalances
    @eachrow! cc.meanbalances begin
      # includes zerosep case by loop design (q = 0)
      :stratum = stratdict[:treatunit]
    end
  end

  cc.stratifier = stratifier

  return cc, nothing # labels should come from elsewhere
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
function variablestrat!(
  cc, dat, var;
  qtes = [0, 0.25, 0.5, 0.75, 1.0], zerosep = false, timevary = nothing
)

  cc.stratifier = var;

  cc.matches[!, :stratum] = Vector{Int}(undef, nrow(cc.matches));
  cc.meanbalances[!, :stratum] = Vector{Int}(undef, nrow(cc.meanbalances));

  c1 = dat[:, cc.treatment] .== 1;

  missingpresent = false;

  # if separate zero
  zs = 0;

  if isnothing(timevary)
    timevary = get(cc.timevary, var, false)
  end

  if !timevary
    udf = unique(@view(dat[c1, :]), [cc.id, var], view = true);
    udict = Dict(udf[!, cc.id] .=> udf[!, var]);
    
    xvec = udf[!, var];
    missingpresent = any(ismissing.(udf[!, var])) # update value
    # if there are missing values in the data, just skip them
    if missingpresent
      xvec = disallowmissing(xvec[.!ismissing.(xvec)]);
    end
    X = sort(quantile(xvec, qtes));
    Xlen = length(X);

    @eachrow! cc.matches begin
      # includes zerosep case by loop design (q = 0)
      dictval = udict[:treatunit]
      # new stratum if missing
      :stratum = !ismissing(dictval) ? assignq(dictval, X, Xlen) : Xlen
    end

    @eachrow! cc.meanbalances begin
      # includes zerosep case by loop design (q = 0)
      dictval = udict[:treatunit]
      # new stratum if missing
      :stratum = !ismissing(dictval) ? assignq(dictval, X, Xlen) : Xlen
    end
    
  elseif timevary # do at time of treatment
    
    X = sort(quantile(@views(dat[c1, var])));
    Xlen = length(X);

    udf = unique(@view(dat[c1, :]), [cc.id, var], view = true);

    # treatment events
    udict = Dict{Tuple{Int64, Int64}, eltype(udf[!, var])}();
    @eachrow udf begin
      udict[($(cc.t), $(cc.id))] = $var
    end
    
    if zerosep
      c2 = dat[:, var] .> 0;
      X = sort(quantile(@views(dat[c2 .& c1, var])));
      zs = 1;
    end

    # matches
    @eachrow! cc.matches begin
      # includes zerosep case by loop design (q = 0)
      :stratum = assignq(udict[(:treattime, :treatunit)], X, Xlen) + zs
    end

    # meanbalances
    @eachrow! cc.meanbalances begin
      # includes zerosep case by loop design (q = 0)
      :stratum = assignq(udict[(:treattime, :treatunit)], X, Xlen) + zs
    end
  end

  stratlabels = label_variablestrat(
    string.(round.(X, digits = 2));
    missingpresent = missingpresent,
    zerosep = zerosep
  )

  return cc, stratlabels
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

"""
    combostrat!(cc, dat, vars::Vector{Symbol}; varslabs = nothing)

Stratify based on the combinations of one or more variables. Strata are formed directly from the variable values.
"""
function combostrat!(cc, dat, vars::Vector{Symbol}; varslabs = nothing)

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
  
  varn = ""
  for (i, var) in enumerate(vars)
    if i < length(vars)
      varn = varn * string(var) * " x "
    else
      varn = varn * string(var)
    end
  end
  
  cc.stratifier = Symbol(varn);
  
  cc.matches[!, :stratum] = Vector{Int}(undef, nrow(cc.matches));
  cc.meanbalances[!, :stratum] = Vector{Int}(undef, nrow(cc.meanbalances));
  
  # matches
  @eachrow! cc.matches begin
    :stratum = udict[(:treattime, :treatunit)]
  end
  
  # meanbalances
  @eachrow! cc.meanbalances begin
    :stratum = udict[(:treattime, :treatunit)]
  end
  
  stratathatexist = unique(cc.meanbalances.stratum);

  stratlabels = Dict{Int, String}();
  for (k, v) in stratmap
    if stratmap[k] ∈ stratathatexist
      stratlabels[v] = combostratlab(vars, k, varslabs)
    end
  end
  
  return cc, stratlabels
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