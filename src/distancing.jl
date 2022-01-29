## distancing.jl

function distances_allocate!(matches, flen, covnum)
  Threads.@threads for i in eachindex(matches)
  # for (i, tob) in enumerate(tobsvec)
    
    tob = @views matches[i];
    
    # at least one f is valid
    valids = vec(sum(tob.mus, dims = 2) .> 0);

    matches[i] = @set tob.distances = [
      fill(Inf, flen, sum(valids)) for _ in 1:covnum + 1
    ];

    # matches[i] = @set tob.distances = Vector{Matrix{Float64}}(
    #   undef, flen, sum(valids)
    # );
    # validmus = permutedims(tob.mus[valids, :]);

    # _distances_allocate!(matches[i].distances, validmus, covnum)
    
    # matches[i] = @set tob.mudistances = MatchDist(undef, sum(valids));
    # fill_mudists!(tobsvec[i].mudistances, emus, efsets, covnum);
  end
  return matches
end

function _distances_allocate!(matches_i_distances, validmus, covnum)
  for (s, e) in enumerate(validmus)
    if e
      matches_i_distances[s] = Vector{Float64}(undef, covnum + 1);
    end
  end
  return matches_i_distances[s]
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

function distances_calculate!(
  matches, observations, ids,
  tg, rg, fmin, Lmin, Lmax, Σinvdict; sliding = false
)
  @inbounds Threads.@threads for i in eachindex(observations)
    ob = observations[i];

    @unpack mus, distances = matches[i];

    valids = vec(sum(mus, dims = 2) .> 0);
    validunits = @views ids[valids];
    validmus = @views mus[valids, :];

    if length(validunits) == 0
      continue
    end

    γcs = eachcol(tg[ob]);
    γrs = eachrow(tg[ob]);
    γtimes = rg[ob];
    mahas = Vector{Float64}(undef, length(γtimes));

    distantiate!(
      distances, mahas,
      eachrow(validmus), validunits,
      ob[1], Σinvdict,
      γcs, γrs, γtimes, tg,
      fmin, Lmin, Lmax;
      sliding = sliding
    );
  end
  return matches
end

# match distances

matchwindow(f, tt, pomin, pomax) = (tt + f) + pomin : (tt + f) + pomax;

function distantiate!(
  distances, mahas,
  validmuscols, validunits,
  tt, Σinvdict,
  γcs, γrs, γtimes, tg,
  fmin, Lmin, Lmax;
  sliding = false
)

  # the sliding method would be for the pretreatment matching period ONLY

  for (m, (unit, muscol)) in enumerate(
    zip(
      validunits, validmuscols
      )
    )
    g = tg[(tt, unit)];

    mahadistancing!(
      mahas, Σinvdict, γrs, eachrow(g), γtimes
    );

    # cnt: since φ will track 1:31, and we will have only those that exist 
    __distantiate!(
      distances, m,
      muscol,
      mahas, tt, γcs, eachcol(g), γtimes, Σinvdict, fmin, Lmin, Lmax;
      sliding = sliding
    )

  end

  return distances
end

function __distantiate!(
  distances, m,
  muscol,
  mahas, tt, γcs, gcs, γtimes, Σinvdict, fmin, Lmin, Lmax;
  sliding = false
)

  # (φ, fb) = collect(enumerate(muscol))[1]

  for (φ, fb) in enumerate(muscol)
    if fb

      # mahalanobis distance
      # each mahalanobis() call is costly, so do calculations in outer look and average (better to preallocate mahas vector...)

      if sliding
        error("unfinished")
        # this should grab the pretreatment crossover window
        fw = matchwindow(φ + fmin - 1, tt, Lmin, Lmax);
      else
        # fixed window that ought to be based on the crossover definition
        # but is selected manually
        fw = Lmin + tt : tt + Lmax;
      end

      distances[1][φ, m] = mahaveraging(mahas, γtimes, fw)

      # caliper distances
      for (c, (γc, gc)) in enumerate(zip(γcs, gcs))
        # this is the match distance for an f, for a covar
        distances[c + 1][φ, m] = caldistancing(
            Σinvdict, γc, gc, γtimes, fw, c
        );
      end
    end
  end

  return distances
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
