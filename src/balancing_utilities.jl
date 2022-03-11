# balancing utilities

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
        Lstd[(l, covar)] = inv(std(
          skipmissing(lvals[!, covar]); corrected = true)
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
        Lstd[(l, covar)] = inv(std(lvals[!, covar]; corrected = true))
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
      meanbalances[!, covar] = Vector{Vector{Vector{Union{Float64, Missing}}}}(undef, length(matches))
    else 
      meanbalances[!, covar] = Vector{Vector{Union{Float64, Missing}}}(undef, length(matches))
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
  for covar in covariates
    if timevary[covar]
      balrw[covar] = [Vector{Union{Missing, Float64}}(missing, Len) for _ in 1:fpresent]
    else
      balrw[covar] = Vector{Union{Missing, Float64}}(missing, fpresent)
    end
  end
  return balrw
end

function getfunion!(funion, matches_i_mus)
  for j in eachindex(funion)
    funion[j] = any(matches_i_mus[:, j])
  end
end