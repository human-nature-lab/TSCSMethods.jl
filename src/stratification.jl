# stratified.jl

#=
all of the below essentially follow from splitting matches

just execute on each group of a grouped dataframe,
and append output
=#

function stratify!(stratfunc::Function, stratargs)

  cc, labels = stratfunc(stratargs...)
  grandbalance!(cc)

  treatednums!(cc::AbstractCICModel)

  return cc, labels
end

"""
assign treatment numbers in each category
"""
function treatednums!(cc::cicmodel)
  cc.treatednum = Dict{Int64, Int64}();
  for s in unique(cc.meanbalances.stratum)
    
    cc.treatednum[s] = nrow(
      unique(@view cc.meanbalances[!, [:treattime, :treatunit]][cc.meanbalances.stratum .== s, :])
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

  missingpresent = false;

  # if separate zero
  zs = 0;

  if isnothing(timevary)
    timevary = get(cc.timevary, var, false)
  end

  if !timevary
    udf = unique(dat, [:fips, var], view = true);
    udict = Dict(udf[!, :fips] .=> udf[!, var]);
    
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

    # meanbalances
    @eachrow! cc.meanbalances begin
    # includes zerosep case by loop design (q = 0)
    dictval = udict[:treatunit]
    # new stratum if missing
    :stratum = !ismissing(dictval) ? assignq(dictval, X, Xlen) : Xlen
    end
    
  elseif timevary # do at time of treatment
    c1 = dat[:, cc.treatment] .== 1;
    
    X = sort(quantile(@views(dat[c1, var])));
    Xlen = length(X);

    udf = unique(@view(dat[c1, :]), [:fips, var], view = true);

    # treatment events
    udict = Dict{Tuple{Int64, Int64}, eltype(udf[!, var])}();
    @eachrow udf begin
      udict[(:running, :fips)] = $var
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
