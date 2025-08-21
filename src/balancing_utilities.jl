# balancing utilities

# Buffer pools for BalanceData optimization
const BALANCE_FLOAT_POOL = Ref{Union{Nothing, Channel{Vector{Float64}}}}(nothing)
const BALANCE_BIT_POOL = Ref{Union{Nothing, Channel{BitVector}}}(nothing)

function __init_balance_pools__()
    if BALANCE_FLOAT_POOL[] === nothing
        pool = Channel{Vector{Float64}}(Threads.nthreads() * 4)
        for _ in 1:(Threads.nthreads() * 4)
            put!(pool, Vector{Float64}(undef, 100))  # Conservative max size
        end
        BALANCE_FLOAT_POOL[] = pool
    end
    
    if BALANCE_BIT_POOL[] === nothing
        pool = Channel{BitVector}(Threads.nthreads() * 4)
        for _ in 1:(Threads.nthreads() * 4)
            put!(pool, BitVector(undef, 100))
        end
        BALANCE_BIT_POOL[] = pool
    end
end

function get_balance_data(size::Int, fill_missing::Bool = true)
    __init_balance_pools__()
    
    if size > 100 || !isready(BALANCE_FLOAT_POOL[]) || !isready(BALANCE_BIT_POOL[])
        # Pool empty or size too large, allocate directly
        return BalanceData(size, fill_missing), false
    end
    
    values = take!(BALANCE_FLOAT_POOL[])
    is_missing = take!(BALANCE_BIT_POOL[])
    
    resize!(values, size)
    resize!(is_missing, size)
    
    if fill_missing
        fill!(is_missing, true)
    end
    
    return BalanceData(values, is_missing), true
end

function return_balance_data(bd::BalanceData, is_pooled::Bool)
    if is_pooled && length(bd) <= 100
        put!(BALANCE_FLOAT_POOL[], bd.values)
        put!(BALANCE_BIT_POOL[], bd.is_missing)
    end
end

"""
gives (1 / (std of treated units)) for each l in the matching period

outputs dict[t-l, covariate] where t-l is l days prior to the treatment

(adjust this for f, in sliding window F-defined match period case, which is the case...)
"""
function std_treated(model::VeryAbstractCICModel, dat::DataFrame)

  trtobs = unique(dat[dat[!, model.treatment] .== 1, [model.t, model.id]])
  sort!(trtobs, [model.t, model.id])
  idx = Int[];
  L = Int[];

  mmin = minimum(model.L)
  fmax = maximum(model.F)

  for to in eachrow(trtobs)

    c1 = dat[!, model.id] .== to[2];
    ct = (dat[!, model.t] .>= (to[1] - 1 + mmin)) .& (dat[!, model.t] .<= to[1] + fmax);

    append!(idx, findall((c1 .& ct)))
    # L relative to tt, not the actual time
    append!(L, dat[c1 .& ct, model.t] .- to[1])
  end

  allvals = @view dat[idx, model.covariates];

  if any([Missing <: eltype(c) for c in eachcol(allvals)])
    Lset = unique(L);
    Lstd = zeros(Float64, length(Lset));
    # Lstd = Dict{Tuple{Int64, Symbol}, Union{Float64, Missing}}()
    Lstd = Dict{Tuple{Int64, Symbol}, Float64}()

    for l in Lset
      lvals = @view allvals[L .== l, :]
      for covar in model.covariates
        Lstd[(l, covar)] = 1.0 / std(
          skipmissing(lvals[!, covar]); corrected = true
        )
      end
    end
  else
    Lset = unique(L);
    Lstd = zeros(Float64, length(Lset));
    Lstd = Dict{Tuple{Int64, Symbol}, Float64}()
    
    for l in Lset
      lvals = @view allvals[L .== l, :]
      for covar in model.covariates
        Lstd[(l, covar)] = 1.0 / std(lvals[!, covar]; corrected = true)
      end
    end
  end
  
  return Lstd 
end

"""
    allocate_meanbalances!(model)

Prepare the meanbalances DataFrame for a model. Meanbalances has number of rows == observations.
"""
function allocate_meanbalances!(model)

  @unpack observations, matches, ids, meanbalances = model;
  @unpack covariates, timevary = model;
  @unpack t, id, treatment = model;
  @unpack F, L = model;

  Len = length(L); Flen = length(F);

  # meanbalances = DataFrame(
  #   treattime = [ob[1] for ob in model.observations],
  #   treatunit = [ob[2] for ob in model.observations],
  #   matchunitsets = [mm for mm in model.matchunits]
  # )

  # need to check observations to see if there are any matches left

  meanbalances[!, :fs] = Vector{Vector{Bool}}(undef, length(matches));

  for covar in covariates
    if timevary[covar]
      meanbalances[!, covar] = Vector{Vector{BalanceData}}(undef, length(matches))
    else 
      meanbalances[!, covar] = Vector{BalanceData}(undef, length(matches))
    end
  end

  _fill_meanbalances!(
    meanbalances, matches, Len, covariates, timevary, Flen
  );

  return model
end

function _fill_meanbalances!(
  meanbalances, matches, Len, covariates, timevary, Flen
)

  # (i, balrw) = collect(enumerate(eachrow(meanbalances)))[1]

  for (i, balrw) in enumerate(eachrow(meanbalances))
    
    # we want to create a vector for f, so long as at least one match unit allows it
    balrw[:fs] = Vector{Bool}(undef, Flen);
    getfunion!(balrw[:fs], matches[i].mus);
    fpresent = sum(balrw[:fs]);
  
    __fill_meanbalances!(
      balrw, fpresent, Len, covariates, timevary
    )
  end
  return meanbalances
end

function __fill_meanbalances!(
  balrw, fpresent, Len, covariates, timevary
)
  pooled_data = []  # Track pooled data for cleanup
  
  try
    for covar in covariates
      if timevary[covar]
        covar_data = BalanceData[]
        for _ in 1:fpresent
          bd, is_pooled = get_balance_data(Len, true)
          push!(covar_data, bd)
          push!(pooled_data, (bd, is_pooled))
        end
        balrw[covar] = covar_data
      else
        bd, is_pooled = get_balance_data(fpresent, true)
        balrw[covar] = bd
        push!(pooled_data, (bd, is_pooled))
      end
    end
  catch e
    # Clean up any allocated data on error
    for (bd, is_pooled) in pooled_data
      return_balance_data(bd, is_pooled)
    end
    rethrow(e)
  end
  
  return balrw
end

function getfunion!(funion, matches_i_mus)
  for j in eachindex(funion)
    funion[j] = any(matches_i_mus[:, j])
  end
end