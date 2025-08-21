# match!.jl

"""
    default_treatmentcategories(x)

(For use as input to match!)
Default treatment history categories.
We look at total count in the pre-treatment crossover
period for both the treated unit and its (potential) match.
If the totals, for each, in the same period fall into the same cateory,
we allow the match.

Here the categories are either not treated, or treated. Ths works for single treatments, and prevents matches from being treated at all during the pretreatment crossover window.
"""
function default_treatmentcategories(x)
  return if x == 0
    0
  else 1
  end
end

"""
    match!(
      model::AbstractCICModel, dat::DataFrame;
      treatcat = nothing, exposure = nothing
    )
  
Perform matching for treatment events, using Mahalanobis distance matching. Additionally, calculate standardized Euclidean distances for the individual covariates are specified.
"""
function match!(
  model::AbstractCICModel, 
  dat::DataFrame;
  treatcat::Function = default_treatmentcategories,
  exposure::Union{Nothing, Symbol} = nothing,
  variancesonly::Bool = true
)::AbstractCICModel

  # using Parameters, Accessors
  # import TSCSMethods:eligibility!,distances_allocate!,samplecovar,distances_calculate!,window_distances!,distaveraging!,rank!
  # treatcat = default_treatmentcategories
  # exposure = nothing
  # variancesonly = true
  # sliding = false

  # Input validation
  if nrow(dat) == 0
    throw(ArgumentError("Input data cannot be empty"))
  end
  
  @unpack observations, matches, ids = model;
  @unpack F, L, id, t, treatment, covariates = model;
  
  # Validate that required columns exist in data
  required_cols = [id, t, treatment]
  append!(required_cols, covariates)
  if !isnothing(exposure)
    push!(required_cols, exposure)
  end
  
  missing_cols = setdiff(required_cols, Symbol.(names(dat)))
  if !isempty(missing_cols)
    throw(ArgumentError("Missing required columns in data: $(missing_cols)"))
  end
  
  if length(observations) == 0
    @warn "No treatment observations found - matching may not be meaningful"
  end
  
  flen = length(F);
  
  # bounds for outcome window
  fmin, fmax = extrema(F)

  # bounds for covariate matching window
  Lmin, Lmax = extrema(L)

  covnum = length(covariates);

  # outcome should be column one
  cdat = Matrix(dat[!, covariates]);

  # eligibility handling
  tg, rg = eligibility!(
    matches, observations, cdat,
    ids, treatcat,
    dat[!, t], dat[!, id], dat[!, treatment],
    fmin, fmax, Lmin; exposure = exposure
  );

  # distance between matches and treated
  distances_allocate!(matches, flen, covnum);

  Σinvdict = samplecovar(dat[!, t], cdat; variancesonly = variancesonly);

  distances_calculate!(
    matches, observations, ids, covariates, tg, rg, fmin, Lmin, Lmax, Σinvdict
  );

  # make distances dimensions match the updated mus, by dropping
  # ineligible units
  # (this only works with fixed window matching)
  for i in eachindex(matches)
    tob = @views matches[i]
    matches[i] = @set tob.distances = tob.distances[.!isinf.(tob.distances[:, 1]), :];
  end

  rank!(matches, flen);

  return model
end

function eligibility!(
  matches, observations, cdat::Matrix{Float64},
  ids, treatcat,
  dat_t, dat_id, dat_trt,
  fmin, fmax, Lmin; exposure = nothing
)
  # eligibility handling
  if isnothing(exposure)
    tg, rg, trtg = make_groupindices(
      dat_t, dat_trt,
      dat_id, ids,
      fmax, Lmin,
      cdat
    );
    
    # tvec = dat[!, t];
    # treatvec = dat[!, treatment];
    # idvec = dat[!, id];
    # uid = ids; fmax, Lmin, cdat;

    eligiblematches!(
      observations, matches,
      rg, trtg, ids, fmin, fmax, treatcat
    );
  else
    tg, rg, trtg, exg = make_groupindices(
      dat_t, dat_trt,
      dat_id, ids,
      fmax, Lmin,
      cdat;
      exvec = exposure
    );

    eligiblematches!(
      observations, matches,
      rg, trtg, ids, fmin, fmax, treatcat, exg = exg
    );
  end

  return tg, rg
end

function eligibility!(
  matches, observations, cdat::Matrix{Union{Missing, Float64}},
  ids, treatcat,
  dat_t, dat_id, dat_trt,
  fmin, fmax, Lmin; exposure = nothing
)

  # dat_t = dat[!, t]
  # dat_id = dat[!, id]
  # dat_trt = dat[!, treatment]

  # eligibility handling
  if isnothing(exposure)
    tg, rg, trtg = make_groupindices(
      dat_t, dat_trt,
      dat_id, ids,
      fmax, Lmin,
      cdat
    );
    
    # tvec = dat[!, t];
    # treatvec = dat[!, treatment];
    # idvec = dat[!, id];
    # uid = ids; fmax, Lmin, cdat;

    eligiblematches!(
      observations, matches,
      tg, rg, trtg, ids, fmin, fmax, treatcat
    );
  else
    tg, rg, trtg, exg = make_groupindices(
      dat_t, dat_trt,
      dat_id, ids,
      fmax, Lmin,
      cdat;
      exvec = exposure
    );

    eligiblematches!(
      observations, matches,
      tg, rg, trtg, ids, fmin, fmax, treatcat, exg = exg
    );
  end

  return tg, rg
end

# Handle the case where there are no covariates (empty covariates array)
function eligibility!(
  matches, observations, cdat::Matrix{Union{}},
  ids, treatcat,
  dat_t, dat_id, dat_trt,
  fmin, fmax, Lmin; exposure = nothing
)
  # When no covariates, create proper Float64 matrix with correct row count but 0 columns
  # The row count should match the data length
  n_obs = length(dat_t)
  empty_cdat = Matrix{Float64}(undef, n_obs, 0)
  return eligibility!(
    matches, observations, empty_cdat,
    ids, treatcat,
    dat_t, dat_id, dat_trt,
    fmin, fmax, Lmin; exposure = exposure
  )
end
