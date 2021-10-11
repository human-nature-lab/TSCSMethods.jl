# matching.jl

"""
give fs for a position
"""
function exbound(pos, mrnge, frnge)
  fs = Vector{Int64}();
  lm = length(mrnge)
  for s in 1:length(frnge)
    if ((pos >= s) & (pos < (s + lm))) # inside bounds? exclude
      append!(fs, s + fmin - 1)
    end
  end
  return fs
end

function rankvect(A)
  vect = zeros(Int64, length(A))
  sp = sortperm(A)
  rnk = 1:length(A)
  begin
    for r in rnk
      vect[rnk[sp][r]] = r
    end
  end
  return vect
end

"""
x is most rapid, then y, then z
"""
subi(y, x, z, Y, X) = (z - 1) * X * Y + X * (y - 1) + x

"""
inverted covariance matrix for mahalanobis distance (all units at t)

inverted sqrt(vars) for balance score calculations (treated units at t)
"""
function calculate_sample_Σs!(
  ut, Σinvdict, σdict, dt, c_treatment, cmat,
  variancesonly::Bool
)
  
  for i = eachindex(ut)
    uti = ut[i]
    c_t = dt .== uti
    Σ = cov(cmat[c_t, :])

    if variancesonly
      Σ[Not(diagind(Σ))] .= 0
    end

    c2 = c_t .& c_treatment;
    if (sum(c2) > 1)
      Σ_treated = cov(cmat[c2, :])
      σdict[uti] = 1 ./ sqrt.(diag(Σ_treated)) # for balance score
    end
    Σinvdict[uti] = pinv(Σ)
  end
  return Σinvdict, σdict
end

function samplecovar(
  dat, uid, covariates, t, id, treatment;
  variancesonly = true
)

  sdat = dat[!, vcat(t, covariates, id, treatment)]

  ut = unique(sdat[:, t])
  cmat = Matrix(sdat[!, covariates])

  c_treatment = sdat[!, id] .∈ Ref(uid)

  Σinvdict = Dict{Int64, Matrix{Float64}}();
  σdict = Dict{Int64, Vector{Float64}}();

  calculate_sample_Σs!(
    ut, Σinvdict, σdict, sdat[!, t], c_treatment, cmat,
    variancesonly
  )
  return Σinvdict, σdict
end

function distance_sums!(
  tmin, Σinvdict, covariates, c1tusdf, c1ousdf, LF, Lrnge, caldists, mdist, lcnt
)
  
  for l = eachindex(Lrnge) # mmin:mmax
    ltime = LF[l]
    if ltime >= tmin
      lcnt[1] += 1

      il = Vector{Float64}(undef, length(covariates))
      jl = similar(il)
      
      fillrow!(
        il, jl, covariates, c1tusdf, c1ousdf, lcnt
      )

      mdist[1] += mahalanobis(il, jl, Σinvdict[LF[l]])
      distance_sums_inner!(caldists, covariates, il, jl, Σinvdict, LF, l)
      
    end
  end
  return lcnt, mdist, caldists
end

function fillrow!(il, jl, covariates, c1tusdf, c1ousdf, lcnt)
  for (cnum, covar) = enumerate(covariates)
    il[cnum] = c1tusdf[lcnt[1], covar] # row for iu at l
    jl[cnum] = c1ousdf[lcnt[1], covar] # row for ju at l

  end
return il,jl
end

function distance_sums_inner!(caldists, covariates, il, jl, Σinvdict, LF, l)
    for cnum = eachindex(covariates)
        caldists[cnum] += weuclidean(
          il[cnum], jl[cnum],
          Σinvdict[LF[l]][cnum,cnum]
        )
    end
    return caldists
end

function matching_barrier_2!(
  possibles,
  tusdf, ousdf,
  iu, ju, it,
  luid, lff, lmm, tmin,
  LF,
  fmin, fmax, Lrnge,
  Σinvdict,
  covariates,
  treatment,
  t,
  distances,
  i, j
)
  # then we simply index LF
  for k = eachindex(fmin:fmax)
    ijk = subi(j, k, i, luid, lff);

    f = (k - 1) + fmin;

    possibles[ijk, :treattime] = it;
    possibles[ijk, :treatunit] = iu;
    possibles[ijk, :matchunit] = ju;
    possibles[ijk, :f] = f;

    # get the pre-treatment period outcome value
    # minimum will be at least the treatment date in data
    # need to exclude if it-1 DNE

    #= if the pre-treatment datapoint does not exist, skip to next k
    (should really just skip to next ju, but here to store info before skipping, this should be a rare case)

    could also store the pre-treatment average here, instead of t - 1

    quicker w/o search, assumes no missing entries, same start date for all units
    =#

    tstrt = max((it + fmin) + mmin, tmin) # begin of data or earliest L
    # this is the first time point in sdf
    # assumes same start date for all units
    trtdx = it - tstrt + 1 # first index is 1

    if trtdx < 2 # treatment must be present, and day before treatment
      # possibles[ijk, :poss] = false; # default already
      continue
    end

    Lf = LF[k : ((k - 1) + lmm)];
    c1j = ousdf[!, t] .∈ Ref(Lf);
    # c2a = sum(c1j) >= min_match_len # minimum match period length
    c2 = !any(e -> e == 1, ousdf[c1j, treatment]);
    
    if c2
      possibles[ijk, :possible] = true;

      # do distance calculations
      if distances

        c1ousdf = @view ousdf[c1j, covariates];
        c1i = tusdf[!, t] .∈ Ref(Lf);
        c1tusdf = @view tusdf[c1i, covariates];
        
        mdist = [0.0]
        caldists = zeros(Float64, length(covariates));
        lcnt = [0]
        
        distance_sums!(
          tmin, Σinvdict,
          covariates, c1tusdf, c1ousdf, LF, Lrnge, caldists, mdist, lcnt
        )
        fillcaldists!(possibles, covariates, caldists, lcnt[1], ijk);
        
        possibles[ijk, :mdist] = mdist[1] / lcnt[1];
        
      end
    end
  end
  return possibles
end

function fillcaldists!(possibles, covariates, caldists, lcnt, ijk)
  for (cnum, calvar) in enumerate(covariates)
    possibles[ijk, calvar] = caldists[cnum] / lcnt
  end
  return possibles
end

function matching_barrier_1!(
  possibles,
  sdf, tusdf, iu, it, uid,
  fmin, fmax, mmin, mmax,
  luid, lff, lmm, tmin,
  Σinvdict, covariates, treatment, t, id,
  distances,
  # tref,
  varset,
  i
)

  for j = eachindex(uid)
    ju = uid[j];
    if ceil(iu, digits = -3) != ceil(ju, digits = -3) # voting-fips-specific condition

      LF = (it + fmin) + mmin : (it + fmax) + mmax;
      ousdf = @view sdf[sdf[!, id] .== ju, varset];

      matching_barrier_2!(
        possibles,
        tusdf, ousdf,
        iu, ju, it,
        luid, lff, lmm, tmin,
        LF,
        fmin, fmax, mmin:mmax,
        Σinvdict,
        covariates,
        treatment,
        # outcome,
        t,
        distances,
        # tref,
        i, j
      );
    end
  end
  
  return possibles
end

# this seems like the best one
function getmatches!(
  possibles,
  rid, rt, uid,
  fmin, fmax, mmin, mmax,
  lmm, luid, lff, tmin,
  Σinvdict,
  covariates,
  dat,
  varset,
  treatment, t, id;
  distances = false,
  # pretreat = -1
)
  
  @inbounds Threads.@threads for i = eachindex(rid)
  # for i = 1:10 #eachindex(rid)
  
    it = rt[i]; # could filter based on something like < rt[i] + fmax
    sdf = @view dat[
      (dat[!, t] .<= it + fmax) .& (dat[!, t] .>= it + fmin + mmin),
      varset,
    ];
    # keep data only from it + fmin + mmin (50 days before first f) up until it + fmax (40 days after treatment event)

    iu = rid[i]; # need tou for distances
    tusdf = @view sdf[sdf[!, id] .== iu, :];
    
    # each thread is overwriting the others' work?
    matching_barrier_1!(
      possibles,
      sdf, tusdf, iu, it, uid,
      fmin, fmax, mmin, mmax,
      luid, lff, lmm, tmin,
      Σinvdict, covariates, treatment, t, id,
      distances,
      # pretreat,
      varset,
      i
    )

  end

  return possibles
end

function processmatches!(cic::cicmodel)
  cic.matches = cic.matches[cic.matches[!, :possible], :];
  select!(cic.matches, Not(:possible))

  sort!(cic.matches, [:treattime, :f, :treatunit, :mdist])

  @linq cic.matches = cic.matches |>
    groupby([:treattime, :f, :treatunit]) |>
    transform!(
      :matchunit => eachindex => :rank, nrow => :numpossible
    )
  return cic
end

function setupmatches(dat, t, id, treatment, mmin, mmax, fmin, fmax)

  # begin setup
  uid = unique(dat[!, id]);
  lmm = length(mmin:mmax);
  luid = length(uid);
  lff = length(fmin:fmax);
  tmin = minimum(dat[!, t]);

  varset = vcat([t, id, treatment, outcome], covariates);

  treatment_points = findall(dat[!, treatment] .== 1);
  rid = @views dat[treatment_points, id];
  rt = @views dat[treatment_points, t];

  L = luid * length(treatment_points) * lff;

  possibles = DataFrame(
    treattime = zeros(Int64, L),
    treatunit = zeros(Int64, L),
    matchunit = zeros(Int64, L),
    f = zeros(Int64, L),
    possible = [false for i in 1:L]
  );

  # add caliper vars
  for calvar in vcat(:mdist, covariates)
    possibles[!, calvar] .= 0.0;
  end

  # get sample variances
  Σinvdict, σd = samplecovar(dat, unique(rid), covariates, t, id, treatment);
  return possibles, uid, lmm, luid, lff, tmin, rid, rt, varset, Σinvdict
end


"""
    match!(cic::cicmodel, dat::DataFrame; distances = true)
  
Perform matching for treatment events, using Mahalanobis distance matching. Optionally specify that standardized Euclidean distances for the individual covariates are specified.
"""
function match!(cc::cicmodel, dat::DataFrame; distances = true)

  cc.matches,
  uid, lmm, luid, lff,
  tmin, rid, rt, varset, Σinvdict = setupmatches(
    dat, cc.t, cc.id, cc.treatment,
    cc.mmin, cc.mmax, cc.fmin, cc.fmax
  );

  getmatches!(
    cc.matches,
    rid, rt, uid,
    cc.fmin, cc.fmax, cc.mmin, cc.mmax,
    lmm, luid, lff, tmin,
    Σinvdict,
    cc.covariates,
    dat,
    varset,
    cc.treatment, cc.t, cc.id;
    distances = distances,
    # pretreat = cc.reference
  );

  processmatches!(cc)

  return cc
end

function refine(cc::AbstractCICModel, refinementnum::Int)
  return @view cc.matches[cc.matches[!, :rank] .<= refinementnum, :]
end
