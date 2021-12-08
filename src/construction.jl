# construction.jl

function makemodel(
  dat, t, id, treatment, outcome,
  covariates, timevary,
  F, L;
  title = "model",
  estimator = "ATT"
)

  observations, ids = observe(dat[!, t], dat[!, id], dat[!, treatment]);
  
  # tobs = make_tobsvec(length(observations), length(ids));
  tobs = make_matches(length(observations), length(ids), length(F));

  return CIC(
    title = title,
    id = id,
    t = t,
    treatment = treatment,
    outcome = outcome,
    covariates = covariates,
    timevary = timevary,
    F = F, L = L,
    observations = observations,
    ids = ids,
    matches = tobs,
    treatednum = length(observations),
    estimator = estimator
  )
end
