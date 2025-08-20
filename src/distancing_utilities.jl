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
  return mdist / accum
end

# Unified distance averaging functions that handle both missing and non-missing data
# efficiently with compile-time specialization

# Helper functions for missing data detection (compile-time optimized)
@inline _has_missing_data(::Float64) = false
@inline _has_missing_data(::Missing) = true
@inline _has_missing_data(::Union{Float64, Missing}) = true  # Runtime check needed

@inline _is_value_missing(::Float64) = false
@inline _is_value_missing(::Missing) = true

# Version with φ and m parameters (for sliding windows)
function distaveraging!(
  distances, dtots::Vector{Vector{T}}, accums, γtimes, fw, φ, m
) where {T}

  # Setup and recycling - type-stable initialization
  if T <: Union{Float64, Missing}
    # For Union types, we don't pre-initialize to avoid missing value issues
    # accums will be zeroed in the loop when we first encounter valid data
    accums_initialized = false
  else
    # For pure Float64, we can safely initialize
    for ι in eachindex(dtots)
      distances[ι][φ, m] = 0.0
      accums[ι] = 0
    end
    accums_initialized = true
  end

  # Main averaging loop
  for (l, τ) in enumerate(γtimes)
    if τ > maximum(fw) # don't bother with the rest
      break
    elseif (τ >= minimum(fw)) 
      # Handle initialization for Union types on first valid data
      if !accums_initialized && T <: Union{Float64, Missing}
        for ι in eachindex(dtots)
          distances[ι][φ, m] = 0.0
          accums[ι] = 0
        end
        accums_initialized = true
      end
      
      # Process data based on type
      if T <: Union{Float64, Missing}
        # Need to check for missing values
        for u in eachindex(dtots)
          val = dtots[u][l]
          if !_is_value_missing(val)
            distances[u][φ, m] += val
            accums[u] += 1
          end
        end
      else
        # Pure Float64 - no missing check needed
        for u in eachindex(dtots)
          distances[u][φ, m] += dtots[u][l]
          accums[u] += 1
        end
      end
    end
  end

  # Finalize averages
  for ι in eachindex(dtots)
    distances[ι][φ, m] = if accums[ι] == 0
      Inf
    else
      distances[ι][φ, m] / accums[ι]
    end
  end
end

# Version without φ and m parameters (for fixed windows)
function distaveraging!(
  drow, dtots::Vector{Vector{T}}, accums, γtimes, fw
) where {T}

  # Setup and recycling - type-stable initialization
  if T <: Union{Float64, Missing}
    # For Union types, we don't pre-initialize
    accums_initialized = false
  else
    # For pure Float64, we can safely initialize
    for ι in eachindex(dtots)
      drow[ι] = 0.0
      accums[ι] = 0
    end
    accums_initialized = true
  end

  # Main averaging loop
  for (l, τ) in enumerate(γtimes)
    if τ > maximum(fw) # don't bother with the rest
      break
    elseif (τ >= minimum(fw))
      # Handle initialization for Union types on first valid data
      if !accums_initialized && T <: Union{Float64, Missing}
        for ι in eachindex(dtots)
          drow[ι] = 0.0
          accums[ι] = 0
        end
        accums_initialized = true
      end
      
      # Process data based on type
      if T <: Union{Float64, Missing}
        # Need to check for missing values
        for u in eachindex(drow)
          val = dtots[u][l]
          if !_is_value_missing(val)
            drow[u] += val
            accums[u] += 1
          end
        end
      else
        # Pure Float64 - no missing check needed
        for u in eachindex(drow)
          drow[u] += dtots[u][l]
          accums[u] += 1
        end
      end
    end
  end

  # Finalize averages
  for ι in eachindex(dtots)
    drow[ι] = if accums[ι] == 0
      Inf
    else
      drow[ι] / accums[ι]
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
    mdist / accum
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
    caldist / accum
  end
end
##
