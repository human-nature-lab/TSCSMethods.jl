# matching.jl

function observe(datt, datid, dattrt)
  obslen = sum(dattrt)
  v = Vector{Tuple{Int, Int}}(undef, obslen);
  cnt = 0
  for (τ, unit, z) in zip(datt, datid, dattrt)
    if z > 0
      cnt += 1
      v[cnt] = (τ, unit)
    end
  end

  return sort(v), unique(datid)
end

function make_matches(obslen, idlen)
  matches = Vector{Tob}(undef, obslen);
  _make_matches!(matches, obslen, idlen);
  return matches
end

function _make_matches!(matches, obslen, idlen)
  for i in 1:obslen
    matches[i] = Tob(
      mus = fill(false, idlen),
      fs = Vector{Vector{Bool}}(undef, idlen),
      ranks = Dict{Int, Vector{Int}}()
    );
  end;
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
    make_groupindices(tvec, treatvec, idvec, uid, fmin, fmax, mmin, cdat)

Get partially overlapping values for treated observations.
"""
function make_groupindices(tvec, treatvec, idvec, uid, fmin, fmax, mmin, cdat)
  #  14.354192 seconds
  # (3.11 M allocations: 283.229 MiB, 14.99% gc time, 7.10% compilation time)
  STypeMat = SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Base.Slice{Base.OneTo{Int64}}}, false};
  STypeVec = SubArray{Int64, 1, Vector{Int64}, Tuple{Vector{Int64}}, false}

  tts = sort(unique(tvec[treatvec .== 1]));
  
  tidx = Dict{Tuple{Int, Int}, STypeMat}();
  ridx = Dict{Tuple{Int, Int}, STypeVec}();
  tridx = Dict{Tuple{Int, Int}, STypeVec}();

  sizehint!(tidx, length(Iterators.product(tts, uid)));
  sizehint!(ridx, length(Iterators.product(tts, uid)));
  sizehint!(tridx, length(Iterators.product(tts, uid)));

  _makegroupindices(
    tidx, ridx, tridx, tts, uid, fmin, fmax, mmin, tvec, idvec, treatvec, cdat
  )
  
  return tidx, ridx, tridx
end

function _makegroupindices(
  tidx, ridx, tridx, tts, uid, fmin, fmax, mmin, tvec, idvec, treatvec, cdat
)
  @floop for (tt, unit) in Iterators.product(tts, uid)
    @init yesrows = Vector{Bool}(undef, length(tvec))
    getyes!(yesrows, tvec, idvec, tt, fmin, fmax, mmin, unit)

    tidx[(tt, unit)] = @views cdat[yesrows, :];
    ridx[(tt, unit)] = @views tvec[yesrows];
    tridx[(tt, unit)] = @views treatvec[yesrows];
  end
  return tidx, ridx, tridx
end

function getyes!(yesrows, tvec, idvec, tt, fmin, fmax, mmin, unit)
  for k in eachindex(tvec)
    yesrows[k] = ((tvec[k] < tt + fmax) & (tvec[k] >= tt + fmin + mmin)) & (idvec[k] == unit)
  end
  return yesrows
end

"""
      getfset!(ftrue, fmin, fmax, pollution, gt, tt)

All fs are allowable if the matchunit is treated no closer than 31 days (before or after) the treatment.
  => mu treatment at 31th day after treatment means first outcome
      window day is after tu's fmax
  => 41st day after mu treatment is beyond its outcome window
      so it is fully elgible to be a match to a unit with an fmin
      on that day. Then the treatment can happen 10 (fmin), days # #    before.
"""
function getfset!(ftrue, fmin, fmax, pollution, gt, tt);

  for (d, τ) in zip(pollution, gt)
    # tt+fmin:tt+fmax, τ+fmin:τ+fmax
    if d > 0
      if τ == tt
        for φ in eachindex(ftrue); ftrue[φ] = false; end
        return ftrue
      elseif abs(tt - τ) >= fmin + fmax + 1
        continue
      else
        if tt - τ > 0
          # consider tt - τ = 1 -> 1 estimatble f (at fmax)
          # true at fmax-(tt-τ)+1-fmin+1:fmax-fmin+1
          for φ in 1:fmax-(tt-τ)+1-fmin
            ftrue[φ] = false
          end
        elseif tt - τ < 0
          for φ in fmax-(tt-τ)+1-fmin+1:fmax-fmin+1
            ftrue[φ] = false
          end
        end
      end
    end
  end
  return ftrue
end

function obsmatches!(tobs, rg, trtg, uid, tt, tu, fmin, fmax, flen, ftrue)

  for (m, mu) in enumerate(uid) # [1]# (j, mu) in enumerate(uid)
    if tu != mu

      # require that panels are balanced
      # will need to handle missingness case eventually
      # data are padded
      
      # reset ftrue to true for all values
      # (preallocated before loops)
      for l in eachindex(ftrue)
        ftrue[l] = true
      end
      
      pollution = trtg[(tt, mu)];

      if sum(pollution) > 0
        getfset!(
          ftrue, fmin, fmax, pollution, rg[(tt, mu)], tt
        )
        if !any(ftrue)
          tobs.mus[m] = false
          # tobs.fs[m] = ftrue
        else
          tobs.mus[m] = true
          tobs.fs[m] = ftrue
        end
      else
        tobs.mus[m] = true
        tobs.fs[m] = ftrue
      end
    end
  end
  return tobs
end

function distance_allocate!(tobsvec, ids, covnum)
  Threads.@threads for i in eachindex(tobsvec)
  # for (i, tob) in enumerate(tobsvec)
    tob = @views tobsvec[i];
    # compression from uid to eligible matches
    assigned = getassigned(tob.mus, tob.fs);
    emus = ids[assigned]; # eligible matches
    efsets = tob.fs[assigned]; # allowable fs for each eligible match
    
    # preallocate for distances
    tobsvec[i] = @set tob.mudistances = MatchDist(undef, sum(assigned));
    fill_mudists!(tobsvec[i].mudistances, emus, efsets, covnum);
  end
  return tobsvec
end

function getassigned(mus, fs)
  assigned = Vector{Bool}(undef, length(mus))
  _getassigned!(assigned, mus, fs)
  return assigned
end

function _getassigned!(assigned, mus, fs)
  for i in eachindex(mus)
    assigned[i] = isassigned(fs, i)
  end
  return assigned
end

# preallocate
function fill_mudists!(mudists, emus, efsets, covnum)
  for em in eachindex(emus)
    # emu = emus[e]
    efs = efsets[em] # logical rep. of fs that exist for a match
    # a 4-len for each f that exists for the match
    mudists[em] = Vector{Vector{Float64}}(undef, sum(efs))
    _fill_mudists!(mudists[em], covnum)
  end
  return mudists
end

function _fill_mudists!(mudists_em, covnum)
  for j in eachindex(mudists_em)
    mudists_em[j] = Vector{Float64}(undef, covnum + 1)
  end
  return mudists_em
end

function _getmatches!(
  observations, tobsvec,
  rg, trtg, ids, fmin, fmax, flen, ftrue
)
  for (ob, tob) in zip(observations, tobsvec)
      (tt, tu) = ob;
      obsmatches!(tob, rg, trtg, ids, tt, tu, fmin, fmax, flen, ftrue);
  end
  return tobsvec
end

"""
    match!(cic::cicmodel, dat::DataFrame)
  
Perform matching for treatment events, using Mahalanobis distance matching. Additionally, calculate standardized Euclidean distances for the individual covariates are specified.
"""
function match!(model::AbstractCICModel, dat)

  @unpack observations, matches, ids = model;
  @unpack F, L, id, t, treatment, covariates = model;
  
  flen = length(F);
  fmin = minimum(F); fmax = maximum(F);
  mmin = minimum(L); mmax = maximum(L);
  covnum = length(covariates);

  cdat = Matrix(dat[!, covariates]);
  
  # model = allocate_matches!(observations, ids)
  #   0.030003 seconds (17.26 k allocations: 68.326 MiB)
  #   2.460215 seconds (17.26 k allocations: 68.326 MiB, 99.39% gc time)

  tg, rg, trtg = make_groupindices(
    dat[!, t], dat[!, treatment],
    dat[!, id], ids,
    fmin, fmax, mmin,
    cdat
  );
  
  GC.gc();

  ftrue = Vector{Bool}(undef, flen);
  _getmatches!(
    observations, matches,
    rg, trtg, ids, fmin, fmax, flen, ftrue
  );

  # distance prealloation
  distance_allocate!(matches, ids, covnum);

  # tobsvec[1].mudistances[1]
  Σinvdict = samplecovar(dat, covariates, t, id, treatment);
  
  GC.gc();

  # 175.818641 seconds (1.79 G allocations: 578.938 GiB, 60.88% gc time)
  _distances!(
    matches, observations, ids, tg, rg, fmin, mmin, mmax, Σinvdict
  );

  rank!(matches, flen)

  return model
end

###
# check allowable fs

# ob = observations[2106]
# tmu = (ob[1], 44007)

# findfirst(ids .== tmu[2])
# findfirst(ids[tobsvec[1].mus] .== tmu[2])

# fset = [true for i in 1:31];
# begin
#   pollution = trtg[tmu];
#   @time getfset!(
#       fset, fmin, fmax, pollution, rg[tmu], ob[1]
#   )
# end
# fset

# ftrue = [true for i in 1:31];
# pollution = trtg[tmu]; gt = rg[tmu]; tt = ob[1];

###

function matchassignments(tobsi, ids; returnefsets = true)
  assigned = getassigned(tobsi.mus, tobsi.fs);
  emus = ids[assigned]; # eligible matches
  
  if returnefsets == false
    return emus
  else
    efsets = tobsi.fs[assigned]; # allowable fs for each eligible match
    return emus, efsets
  end
end

function _distances!(
  tobs, observations, ids,
  tg, rg, fmin, mmin, mmax, Σinvdict
)
  @inbounds Threads.@threads for i in eachindex(observations)
    ob = observations[i];

    emus, efsets = matchassignments(tobs[i], ids)

    if length(emus) == 0
      continue
    end

    γcs = eachcol(tg[ob]);
    γrs = eachrow(tg[ob]);
    γtimes = rg[ob];
    mahas = Vector{Float64}(undef, length(γtimes));
    distantiate!(
      tobs[i].mudistances, mahas,
      ob[1], Σinvdict,
      γcs, γrs, γtimes, tg,
      emus, efsets,
      fmin, mmin, mmax
    );
  end
  return tobs
end

# match distances

matchwindow(f, tt, mmin, mmax) = (tt + f) + mmin : (tt + f) + mmax;

function distantiate!(
  mudists, mahas,
  tt, Σinvdict,
  γcs, γrs, γtimes, tg,
  emus, efsets,
  fmin, mmin, mmax
)
  for (em, emu) in enumerate(emus)
    # emu = emus[e];
    efs = efsets[em]; # logical rep. of fs that exist for a match
    
    g = tg[(tt, emu)];

    mahadistancing!(
        mahas, Σinvdict, γrs, eachrow(g), γtimes
    );

    # cnt: since φ will track 1:31, and we will have only those that exist 
    __distantiate!(
      mudists[em],
      efs, mahas, tt, γcs, eachcol(g), γtimes, Σinvdict, fmin, mmin, mmax
    )
  end
  return mudists
end

function __distantiate!(
  mudists_em, efs, mahas, tt, γcs, gcs, γtimes, Σinvdict, fmin, mmin, mmax
)
  cnt = 0
  for (φ, fb) in enumerate(efs) # eachindex(mudists[em])
    if fb
      cnt += 1

      fw = matchwindow(φ + fmin - 1, tt, mmin, mmax);
      
      # mahalanobis distance
      # each mahalanobis() call is costly, so do calculations in outer look and average (better to preallocate mahas vector...)
      # mudists[em][cnt][1] = @time mahadistancing(
      #     Σinvdict, eachrow(γ), eachrow(g), γtimes, fw
      # );

      mudists_em[cnt][1] = mahaveraging(mahas, γtimes, fw)

      # caliper distances
      for (c, (γc, gc)) in enumerate(zip(γcs, gcs))
        # this is the match distance for an f, for a covar
        mudists_em[cnt][c+1] = caldistancing(
            Σinvdict, γc, gc, γtimes, fw, c
        );
      end

      # figure out distancing here
      # γ; g

      # nightmare b/c of type instability
      # @time mahadistancing(Σinvdict, xrows, yrows, T, fw);

    end
  end
  return mudists_em
end

function mahadistancing(Σinvdict, xrows, yrows, T, fw)
  mahadist = 0.0
  for (xr, yr, τ) in zip(xrows, yrows, T)
    if τ ∈ fw
      Σ = get(Σinvdict, τ, nothing);
      mahadist += sqrt((xr - yr)' * Σ * (xr - yr));
      # function from Distances.jl is really slow?
      # mahadist += @time mahalanobis(xr, yr, Σ);
    end
  end
  return mahadist / length(T)
end

function mahaveraging(mahas, T, fw)
  mdist = 0.0
  for (m, τ) in zip(mahas, T)
    if τ ∈ fw
      mdist += m
    end
  end
  return mdist / length(mahas)
end

function mahadistancing!(mahas, Σinvdict, xrows, yrows, T)
  for (xr, yr, τ, i) in zip(xrows, yrows, T, 1:length(T))
    Σ = get(Σinvdict, τ, nothing);
    mahas[i] = sqrt((xr - yr)' * Σ * (xr - yr));
    # function from Distances.jl is really slow?
    # mahadist += @time mahalanobis(xr, yr, Σ);
  end
  return mahas
end

function caldistancing(Σinvdict, X, Y, T, fw, c)
  caldist = 0.0
  for (x, y, τ) in zip(X, Y, T)
    if τ ∈ fw
      Σ = get(Σinvdict, τ, nothing);
      caldist += weuclidean(x, y, Σ[c,c])
    end
  end
  return caldist / length(X)
end

function rank!(matches, flen)
  Threads.@threads for i in 1:length(matches)
    tob = @views matches[i];
    _rankmatches!(tob, flen)
  end
  return matches
end

function _rankmatches!(tob, flen)
  # murank: f => mu[mu] order ranking
  # 
  @unpack mus, fs, mudistances, ranks = tob

  positions = Vector{Int}(undef, sum(mus))
  cntp = 0
  for (p, mu) in enumerate(mus)
    if mu
      cntp += 1
      positions[cntp] = p
    end
  end

  mamat = fill(Inf, flen, length(positions));

  _rankmatches!(mamat, positions, fs, mudistances, 1:flen);

  for (φ, ec) in enumerate(eachrow(mamat))
    ranks[φ] = positions[sortperm(ec)] # this will put impossible, Inf, values at end -- so beware -- they will be ranked, rely on other stuff to prevent them from ever mattering
  end
  return tob
end

function _rankmatches!(mamat, positions, fs, mudistances, Φ)
  mcnt = 0
  for m in eachindex(positions)
    mcnt += 1

    fsetm = @views fs[positions][m];
    
    # these are the φs included
    # φs = [eachindex(F)...][fsetm];

    φs = Vector{Int}(undef, sum(fsetm))
    cnt = 0
    for (fb, ff) in zip(fsetm, Φ)
      if fb
        cnt += 1
        φs[cnt] = ff
      end
    end
    
    for (j, φ) in enumerate(φs)
      mamat[φ, m] = mudistances[m][j][1] # index = 1 for maha distance
    end

    # _innerassign!(Q, fpos, mudistances)

  end
  return mamat
end