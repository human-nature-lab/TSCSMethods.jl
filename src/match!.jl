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
  model::AbstractCICModel, dat;
  treatcat::Function = default_treatmentcategories,
  exposure = nothing,
  variancesonly = true
)

  # using Parameters
  # import TSCSMethods:eligibility!,distances_allocate!,samplecovar,distances_calculate!,window_distances!,distaveraging!,rank!
  # treatcat = default_treatmentcategories
  # exposure = nothing
  # variancesonly = true

  @unpack observations, matches, ids = model;
  @unpack F, L, id, t, treatment, covariates = model;
  
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
  )

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

  # remove treated observations with no valid mus
  anymatches = fill(true, length(observations));
  for (i, e) in enumerate(matches)
    anymatches[i] = any(e.mus)
  end

  model = @set model.observations = observations[anymatches];
  model = @set model.matches = matches[anymatches];
  model = @set model.treatednum = length(model.observations);

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
