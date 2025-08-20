## distancing.jl

# Thread-local storage for distance calculation working arrays
# This eliminates repeated allocations in the hot path
mutable struct ThreadLocalDistanceStorage
    dtots_float::Vector{Vector{Float64}}
    dtots_missing::Vector{Vector{Union{Float64, Missing}}}
    accums::Vector{Int}
    max_times::Int
    max_covariates::Int
end

# Global thread-local storage
const THREAD_DISTANCE_STORAGE = Vector{Union{Nothing, ThreadLocalDistanceStorage}}(nothing, Threads.nthreads())

function get_thread_storage(n_times::Int, n_covariates::Int)::ThreadLocalDistanceStorage
    thread_id = Threads.threadid()
    storage = THREAD_DISTANCE_STORAGE[thread_id]
    
    # Initialize or resize storage if needed
    if storage === nothing || storage.max_times < n_times || storage.max_covariates < n_covariates
        new_max_times = max(n_times, storage === nothing ? 0 : storage.max_times)
        new_max_covariates = max(n_covariates, storage === nothing ? 0 : storage.max_covariates)
        
        storage = ThreadLocalDistanceStorage(
            [Vector{Float64}(undef, new_max_times) for _ in 1:(new_max_covariates+1)],
            [Vector{Union{Float64, Missing}}(undef, new_max_times) for _ in 1:(new_max_covariates+1)],
            Vector{Int}(undef, new_max_covariates+1),
            new_max_times,
            new_max_covariates
        )
        THREAD_DISTANCE_STORAGE[thread_id] = storage
    end
    
    return storage
end

function distances_allocate!(matches, flen, covnum; sliding = false)
  Threads.@threads :greedy for i in eachindex(matches)
  # for (i, tob) in enumerate(tobsvec)
    
    tob = @views matches[i];
    
    # at least one f is valid
    valids = vec(sum(tob.mus, dims = 2) .> 0);

    if sliding
      matches[i] = @set tob.distances = [
        fill(Inf, flen, sum(valids)) for _ in 1:covnum + 1
      ];
    else
      matches[i] = @set tob.distances = fill(Inf, sum(valids), covnum+1)
    end
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
  matches, observations, ids, covariates,
  tg, rg, fmin, Lmin, Lmax, Σinvdict; sliding = false
)

  @inbounds Threads.@threads :greedy for i in eachindex(observations)
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
    γrs = eachrow(tg[ob]);
    γtimes = rg[ob];
    
    # Get thread-local pre-allocated storage
    storage = get_thread_storage(length(γtimes), length(covariates))
    
    # Use pre-allocated arrays (resize views if needed)
    has_missing = Missing <: eltype(tg[ob])
    dtots = if has_missing
      # Use missing-compatible arrays, resize to current needs
      for h in 1:(length(covariates)+1)
        if length(storage.dtots_missing[h]) < length(γtimes)
          resize!(storage.dtots_missing[h], length(γtimes))
        end
        fill!(view(storage.dtots_missing[h], 1:length(γtimes)), missing)
      end
      storage.dtots_missing
    else
      # Use Float64 arrays
      for h in 1:(length(covariates)+1)
        if length(storage.dtots_float[h]) < length(γtimes)
          resize!(storage.dtots_float[h], length(γtimes))
        end
      end
      storage.dtots_float
    end
    
    # Use pre-allocated accumulator array
    accums = storage.accums
    ##

    window_distances!(
      distances,
      dtots, accums,
      eachrow(validmus), validunits,
      ob[1], Σinvdict,
      γrs, γtimes, tg,
      fmin, Lmin, Lmax;
      sliding = sliding
    );

    # (this will only work, as-is, for non-sliding window)
    # get the matches for which mahalanobis distance cannot be calculated
    # -- due to missingness.

    @views(mus[valids, :][isinf.(distances[:, 1]), :]) .= false;
    # @views(mus[valids, :][.!isinf.(distances[:, 1]), :])
    # @set matches[i].distances = distances[.!isinf.(distances[:, 1]), :];
    # @reset matches[i].distances = distances[.!isinf.(distances[:, 1]), :];
    
  end
  return matches
end

# match distances

matchwindow(f, tt, pomin, pomax) = (tt + f) + pomin : (tt + f) + pomax;

"""
Assign distances for each window. This could be massively
simplified for a non-sliding window.
"""
function window_distances!(
  distances, dtots, accums,
  validmuscols, validunits,
  tt, Σinvdict,
  γrs, γtimes, tg,
  fmin, Lmin, Lmax;
  sliding = false
)

  # the sliding method would be for the pretreatment matching period ONLY
  # cc = 0

  # validmuscols = eachrow(validmus)
  # (m, (unit, muscol)) = collect(enumerate(
  #   zip(
  #     validunits, validmuscols
  #     )
  #   ))[6]

  for (m, (unit, muscol)) in enumerate(
    zip(
      validunits, validmuscols
      )
    )

    # cc += 1
    g = tg[(tt, unit)];

    alldistances!(dtots, Σinvdict, γrs, eachrow(g), γtimes);

    if sliding
      # MISSING DONE

      # @time mahadistancing!(
      #   mahas, Σinvdict, γrs, eachrow(g), γtimes
      # );

      # cnt: since φ will track 1:31, and we will have only those that exist 
      _window_distances!(
        distances, m,
        muscol,
        dtots, accums,
        tt, γtimes, fmin, Lmin, Lmax;
        sliding = sliding
      )

    else
      # NOT SLIDING
      _window_distances!(
        distances, m,
        muscol,
        dtots, accums,
        tt, γtimes, fmin, Lmin, Lmax;
        sliding = sliding
      );
    end
  end

  return distances
end


"""
N.B. if any of the input xr or yr (the covariate values for the treated and match at some day in the covariate matching window) are missing, then the 
maha distance is missing.
"""
function alldistances!(dtotals, Σinvdict, xrows, yrows, γtimes)
  for (xr, yr, τ, i) in zip(xrows, yrows, γtimes, 1:length(γtimes))
    Σ = get(Σinvdict, τ, nothing);

    # mahalanobis distance
    dtotals[1][i] = sqrt((xr - yr)' * Σ * (xr - yr));

    # individuals covariate distances (for calipers)
    for (j, (xj, yj)) in enumerate(zip(xr, yr))
      if !ismissing(xj) & !ismissing(yj)
        dtotals[j+1][i] = weuclidean(xj, yj, Σ[j, j])
      else
        dtotals[j+1][i] = missing
      end
    end

  end
  return dtotals
end

function _window_distances!(
  distances, m,
  muscol,
  dtots, accums,
  tt, γtimes, fmin, Lmin, Lmax;
  sliding = false
)

  if !sliding
    fw = Lmin + tt : tt + Lmax;
    drow = @views distances[m, :];
    ## recycling
    drow .= 0.0
    accums .= 0
    ##
    distaveraging!(drow, dtots, accums, γtimes, fw);

  else
    error("optioned method is unfinished")
    for (φ, fb) in enumerate(muscol)
      if fb # if the specific f (for the given match) is valid

        # mahalanobis distance
        # each mahalanobis() call is costly, so do calculations in outer look and average (better to preallocate mahas vector...)

        # an alternative would be to use looping to find the indices
        # and then just use mean with skipmissing on the appropriate
        # portion of mahas...
        
        fw = matchwindow(φ + fmin - 1, tt, Lmin, Lmax);
        # fw = if sliding
        #   error("unfinished")
        #   # this should grab the pretreatment crossover window
        #   matchwindow(φ + fmin - 1, tt, Lmin, Lmax);
        # else
        #   # This is a gixed window that ought to be based
        #   # on the crossover definition.
        #   # It is selected manually.
        #   Lmin + tt : tt + Lmax;
        # end
        

        distaveraging!(distances, dtots, accums, γtimes, fw, φ, m);
      end

    end
  end

  return distances
end
