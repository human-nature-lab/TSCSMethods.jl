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
    match!(cic::cicmodel, dat::DataFrame; treatcat = nothing)
  
Perform matching for treatment events, using Mahalanobis distance matching. Additionally, calculate standardized Euclidean distances for the individual covariates are specified.
"""
function match!(
  model::AbstractCICModel, dat;
  treatcat::Function = default_treatmentcategories
)

  @unpack observations, matches, ids = model;
  @unpack F, L, id, t, treatment, covariates = model;
  
  flen = length(F);
  fmin = minimum(F); fmax = maximum(F);
  mmin = minimum(L); mmax = maximum(L);
  covnum = length(covariates);

  cdat = Matrix(dat[!, covariates]);

  tg, rg, trtg = make_groupindices(
    dat[!, t], dat[!, treatment],
    dat[!, id], ids,
    fmin, fmax, mmin,
    cdat
  );

  GC.gc();

  eligiblematches!(
    observations, matches,
    rg, trtg, ids, fmin, fmax, treatcat
  );

  # "eligibility" for units
  # eligibility = similar(matches[1].mus) .* 0;
  # over all treated units
  # => num. times a particular unit is eligible to be a match, for each F
  # in the outcome window
  # for objet in matches
  #   eligibility += objet.mus
  # end

  distances_allocate!(matches, flen, covnum);

  Σinvdict = samplecovar(dat, covariates, t, id, treatment);
  
  GC.gc();

  # 193.266173 seconds
  # (2.17 G allocations: 677.980 GiB, 59.13% gc time, 0.01% compilation time)
  distances_calculate!(
    matches, observations, ids, tg, rg, fmin, mmin, mmax, Σinvdict
  );

  rank!(matches, flen);

  return model
end

function samplecovar(
  dat, covariates, t, id, treatment;
  variancesonly = true
)

  sdat = dat[!, vcat(t, covariates, id, treatment)]

  ut = unique(sdat[:, t])
  cmat = Matrix(sdat[!, covariates])

  # c_treatment = sdat[!, id] .∈ Ref(uid)

  Σinvdict = Dict{Int64, Matrix{Float64}}();
  # σdict = Dict{Int64, Vector{Float64}}();

  calculate_sample_Σs!(
    ut, Σinvdict, sdat[!, t], cmat,
    variancesonly
  )
  return Σinvdict
end

"""
inverted covariance matrix for mahalanobis distance (all units at t)

inverted sqrt(vars) for balance score calculations (treated units at t)
"""
function calculate_sample_Σs!(
  ut, Σinvdict, dt, cmat,
  variancesonly::Bool
)
  
  for i = eachindex(ut)
    uti = ut[i]
    c_t = dt .== uti
    Σ = cov(cmat[c_t, :])

    if variancesonly
      Σ[Not(diagind(Σ))] .= 0
    end

    # c2 = c_t .& c_treatment;
    # if (sum(c2) > 1)
    #   Σ_treated = cov(cmat[c2, :])
    #   σdict[uti] = 1 ./ sqrt.(diag(Σ_treated)) # for balance score
    # end
    Σinvdict[uti] = pinv(Σ)
  end
  return Σinvdict
end
