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
    samplecovar_missingness(cdat)
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

function samplecovar_missingness(cdat)
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

function mahaveraging(mahas::Vector{Union{Float64, Missing}}, T, fw)
  mdist = 0.0
  accum = 0
  for (m, τ) in zip(mahas, T)
    if (τ ∈ fw) & !ismissing(m)
      accum += 1
      mdist += m
    end
  end
  return mdist * inv(accum)
end

function caldistancing(Σinvdict, X, Y, T, fw, c)
  caldist = 0.0
  accum = 1
  for (x, y, τ) in zip(X, Y, T)
    if (τ ∈ fw) & (!ismissing(x)) & (!ismissing(y))
      accum += 1
      Σ = get(Σinvdict, τ, nothing);
      caldist += weuclidean(x, y, Σ[c, c])
    end
  end
  return caldist * inv(accum)
end
##
