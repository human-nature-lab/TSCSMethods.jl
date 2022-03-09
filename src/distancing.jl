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

  end
  return matches
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

function distances_calculate!(
  matches, observations, ids,
  tg, rg, fmin, Lmin, Lmax, Σinvdict; sliding = false
)

  @inbounds Threads.@threads for i in eachindex(observations)
    ob = observations[i];

    @unpack mus, distances = matches[i];

    ## check the set of valid matches for ob
    valids = vec(sum(mus, dims = 2) .> 0);
    validunits = @views ids[valids];
    validmus = @views mus[valids, :];

    # skip to next ob if there are no valid matches
    if length(validunits) == 0
      continue
    end
    ##

    ## setup up for distance calculations
    γcs = eachcol(tg[ob]);
    γrs = eachrow(tg[ob]);
    γtimes = rg[ob];
    
    mahas = if Missing <: eltype(tg[ob])
      Vector{Union{Float64, Missing}}(missing, length(γtimes));
    else
      Vector{Float64}(undef, length(γtimes));
    end
    ##

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

  # cc = 0

  # (m, (unit, muscol)) = collect(enumerate(
  #   zip(
  #     validunits, validmuscols
  #     )
  #   ))[1]

  for (m, (unit, muscol)) in enumerate(
    zip(
      validunits, validmuscols
      )
    )
    # cc += 1
    g = tg[(tt, unit)];

    # MISSING DONE
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

"""
N.B. if any of the input xr or yr (the covariate values for the treated and match at some day in the covariate matching window) are missing, then the 
maha distance is missing.
"""
function mahadistancing!(mahas, Σinvdict, xrows, yrows, T)
  for (xr, yr, τ, i) in zip(xrows, yrows, T, 1:length(T))
    Σ = get(Σinvdict, τ, nothing);

    mahas[i] = sqrt((xr - yr)' * Σ * (xr - yr));
    # function from Distances.jl is really slow?
    # mahadist += @time mahalanobis(xr, yr, Σ);
  end
  return mahas
end

function __distantiate!(
  distances, m,
  muscol,
  mahas, tt, γcs, gcs, γtimes, Σinvdict, fmin, Lmin, Lmax;
  sliding = false
)

  # (φ, fb) = collect(enumerate(muscol))[1]

  cc = 0
  for (φ, fb) in enumerate(muscol)
    cc += 1
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

# function mahadistancing(Σinvdict, xrows, yrows, T, fw)
#   mahadist = 0.0
#   accum = 0
#   for (xr, yr, τ) in zip(xrows, yrows, T)
#     if τ ∈ fw
#       accum += 1
#       Σ = get(Σinvdict, τ, nothing);
#       mahadist += sqrt((xr - yr)' * Σ * (xr - yr));
#       # function from Distances.jl is really slow?
#       # mahadist += @time mahalanobis(xr, yr, Σ);
#     end
#   end
#   return mahadist * inv(accum)
# end
