# distancing_utilities.jl

## covariance calculations for mahalanobis distance matching

function samplecovar(
  dat_t, cdat;
  variancesonly = true
)

  ## setup
  ut = sort(unique(dat_t));
  Σinvdict = Dict{Int64, Matrix{Float64}}();
  ##

  ## handle missingness
  cdat2, dat_t2 = if Missing <: eltype(typeof(cdat))
    samplecovar_missingness(cdat, dat_t)
  else
    cdat, dat_t
  end
  ##
  
  calculate_sample_Σs!(
    ut, Σinvdict, dat_t2, cdat2,
    variancesonly
  )
  return Σinvdict
end

function samplecovar_missingness(cdat, dat_t)
  nomis = fill(true, size(cdat)[1]);
  for (k, r) in enumerate(eachrow(cdat))
    if ismissing(sum(r))
      nomis[k] = false
    end
  end

  # use subset of data where each row has no missing values
  return @views(cdat[nomis, :]), @views(dat_t[nomis])
end

"""
inverted covariance matrix for mahalanobis distance (all units at t)

inverted sqrt(vars) for balance score calculations (treated units at t)
"""
function calculate_sample_Σs!(
  ut, Σinvdict, dat_t, cdat,
  variancesonly::Bool
)
  
  for i = eachindex(ut)
    uti = ut[i]
    c_t = dat_t .== uti;
    Σ = cov(cdat[c_t, :]);

    if variancesonly
      Σ[Not(diagind(Σ))] .= 0
    end

    Σinvdict[uti] = pinv(Σ)
  end
  return Σinvdict
end
##

## distance calculations

function mahaveraging(mahas::Vector{Float64}, T, fw)
  mdist = 0.0
  accum = 0
  for (m, τ) in zip(mahas, T)
    if τ ∈ fw
      accum += 1
      mdist += m
    end
  end
  return mdist * inv(accum)
end

function distaveraging!(
  distances, dtots::Vector{Vector{Union{Float64, Missing}}},
  accums, γtimes, fw, φ, m
)

  for (l, τ) in enumerate(γtimes)
    if τ > maximum(fw) # don't bother with the rest
      break
    elseif (τ >= minimum(fw)) & !ismissing(m) # if at or above bottom (above ruled out already)
      # accum += 1
      for u in eachindex(dists)
        if !ismissing(dtots[u][l])
          distances[u][φ, m] += dtots[u][l]
          accums[u] += 1
        end
      end
    end
  end

  for ι in eachindex(dtots)
    distances[ι][φ, m] = if accums[ι] == 0
      Inf
    else
      distances[ι][φ, m] * inv(accums[ι])
    end
  end
end


function distaveraging!(
  drow, dtots::Vector{Vector{Union{Float64, Missing}}},
  accums, γtimes, fw
)

  for (l, τ) in enumerate(γtimes)
    if τ > maximum(fw) # don't bother with the rest
      break
    elseif (τ >= minimum(fw)) & !ismissing(m) # if at or above bottom (above ruled out already)
      # accum += 1
      for u in eachindex(drow)
        if !ismissing(dtots[u][l])
          drow[u] += dtots[u][l]
          accums[u] += 1
        end
      end
    end
  end

  for ι in eachindex(dtots)
    drow[ι] = if accums[ι] == 0
      Inf
    else
      drow[ι] * inv(accums[ι])
    end
  end
end


function distaveraging!(
  distances, dtots::Vector{Vector{Float64}}, accums, γtimes, fw, φ, m
)

  ## setup and recycling
  for ι in eachindex(dtots)
    distances[ι][φ, m] = 0.0
    accums[ι] = 0
  end
  ##

  for (l, τ) in enumerate(γtimes)
    if τ > maximum(fw) # don't bother with the rest
      break
    elseif (τ >= minimum(fw)) # if at or above bottom (above ruled out already)
      # accum += 1
      for u in eachindex(dists)
        distances[u][φ, m] += dtots[u][l]
        accums[u] += 1
      end
    end
  end

  for ι in eachindex(dtots)
    distances[ι][φ, m] = if accums[ι] == 0
      Inf
    else
      distances[ι][φ, m] * inv(accums[ι])
    end
  end
end

function distaveraging!(
  drow,
  dtots::Vector{Vector{Float64}}, accums, γtimes, fw
)

  ## setup and recycling
  for ι in eachindex(dtots)
    drow[ι] = 0.0
    accums[ι] = 0
  end
  ##

  for (l, τ) in enumerate(γtimes)
    if τ > maximum(fw) # don't bother with the rest
      break
    elseif (τ >= minimum(fw)) # if at or above bottom (above ruled out already)
      # accum += 1
      for u in eachindex(drow)
        drow[u] += dtots[u][l]
        accums[u] += 1
      end
    end
  end

  for ι in eachindex(dtots)
    drow[ι] = if accums[ι] == 0
      Inf
    else
      drow[ι] * inv(accums[ι])
    end
  end
end

function mahaveraging(mahas::Vector{Union{Float64, Missing}}, T, fw)
  mdist = 0.0
  accum = 0
  for (m, τ) in zip(mahas, T)
    if τ > maximum(fw) # don't bother with the rest
      break
    elseif (τ >= minimum(fw)) & !ismissing(m) # if at or above bottom (above ruled out already)
      accum += 1
      mdist += m
    end
  end

  return if accum == 0
    Inf
  else
    mdist * inv(accum)
  end
end

function caldistancing(Σinvdict, X, Y, T, fw, c)
  caldist = 0.0
  accum = 0
  for (x, y, τ) in zip(X, Y, T)
    if (τ > maximum(fw)) & (!ismissing(x)) & (!ismissing(y))
      accum += 1
      Σ = get(Σinvdict, τ, nothing);
      caldist += weuclidean(x, y, Σ[c, c])
    end
  end

  return if accum == 0
    Inf
  else
    caldist * inv(accum)
  end
end
##
