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
  qtes = [0, 0.25, 0.5, 0.75, 1.0], zerosep = false
)

  cc.stratifier = var;

  cc.matches[!, :stratum] = Vector{Int}(undef, nrow(cc.matches));
  cc.meanbalances[!, :stratum] = Vector{Int}(undef, nrow(cc.meanbalances));

  if !cc.timevary[var]
    udf = unique(dat, [:fips, var], view = true);
    
    udict = Dict(udf[!, :fips] .=> udf[!, var]);
    
    X = sort(quantile(udf[!, var], qtes));
    Xlen = length(X);

    @eachrow! cc.matches begin
      # includes zerosep case by loop design (q = 0)
      :stratum = assignq(udict[:treatunit], X, Xlen)
    end

    # meanbalances
    @eachrow! cc.meanbalances begin
      # includes zerosep case by loop design (q = 0)
      :stratum = assignq(udict[:treatunit], X, Xlen)
    end
    
  elseif cc.timevary # do at time of treatment
    c1 = dat[:, cc.treatment] .== 1;
    X = sort(quantile(@views(dat[c1, var])));

    udf = unique(@view(dat[c1, :]), [:fips, var], view = true);

    # treatment events
    udict = Dict{Tuple{Int64, Int64}, Int64}();
    @eachrow udf begin
      udict[(:running, :fips)] = $var
    end
    
    zs = 0;
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

  return cc, label_variablestrat(string.(round.(X, digits = 2)))
end

function label_variablestrat(quantnames)
  labels = Dict{Int, String}();
  sizehint!(labels, length(quantnames)-1);
  for i in 1:length(quantnames) - 1
    labels[i] = quantnames[i] * " to " * quantnames[i+1]
  end
  return labels
end