# boostrapping.jl

# Thread-safe bootstrap matrix pool using channels
const MATRIX_POOL = let
    max_flen, max_iter = 200, 5000
    pool_size = Threads.nthreads() * 2  # Allow some extra capacity
    
    pool = Channel{Tuple{Matrix{Float64}, Matrix{Float64}}}(pool_size)
    
    # Pre-populate pool with matrix pairs
    for _ in 1:pool_size
        put!(pool, (
            Matrix{Float64}(undef, max_flen, max_iter),
            Matrix{Float64}(undef, max_flen, max_iter)
        ))
    end
    
    pool
end

function setup_bootstrap(Flen, iterations)
    # Conservative limits for pooling
    max_flen, max_iter = 200, 5000
    
    if Flen > max_flen || iterations > max_iter
        # Dimensions too large, allocate directly
        boots = zeros(Float64, Flen, iterations)
        tcountmat = zeros(Float64, Flen, iterations)
        return boots, tcountmat
    end
    
    # Try to get matrices from pool (non-blocking)
    if isready(MATRIX_POOL)
        boots_buffer, tcounts_buffer = take!(MATRIX_POOL)
        
        # Return matrices to pool when done (this happens automatically when they go out of scope)
        finalizer(_ -> put!(MATRIX_POOL, (boots_buffer, tcounts_buffer)), boots_buffer)
        
        # Create views of the appropriate size
        boots = @view boots_buffer[1:Flen, 1:iterations]
        tcounts = @view tcounts_buffer[1:Flen, 1:iterations]
        
        # Zero out working areas
        fill!(boots, 0.0)
        fill!(tcounts, 0.0)
        
        return boots, tcounts
    else
        # Pool empty, allocate directly
        boots = zeros(Float64, Flen, iterations)
        tcountmat = zeros(Float64, Flen, iterations)
        return boots, tcountmat
    end
end

"""
        bootstrap!(boots, tcountmat, fblocks, ids, treatdex, iterations)

Perform bootstrapping of the att for each f in the outcome window.
Does so for the given set of treated units and matches contained
in the fblocks.
"""
function bootstrap!(boots, tcountmat, fblocks, ids, treatdex, iterations)
    @inbounds Threads.@threads :greedy for b in 1:iterations
        bootcol = @views boots[:, b]
        tcountcol = @views tcountmat[:, b]
        bootatt!(bootcol, tcountcol, fblocks, ids, treatdex);
    end
    return boots
end

function bootatt!(atts, tcounts, fblocks, ids, treatdex)
    sampcount = getsample(ids, treatdex) # randomness
    _boot!(atts, tcounts, fblocks, sampcount)
end

function _boot!(atts, tcounts, fblocks, sampcount)
    for φ in 1:length(fblocks)
        @unpack matchunits, weightedoutcomes,
        weightedrefoutcomes, treatment = fblocks[φ]

        __boot!(
            atts, tcounts, φ,
            matchunits, weightedoutcomes,
            weightedrefoutcomes, treatment,
            sampcount
        )

        atts[φ] = atts[φ] / tcounts[φ]
    end
end

function __boot!(
    atts, tcounts, φ,
    matchunits, weightedoutcomes,
    weightedrefoutcomes, treatments,
    sampcount
)
    for (munit, wo, wref, trted) in zip(
        matchunits, weightedoutcomes,
        weightedrefoutcomes, treatments
    )
        atts[φ] += (wo + wref) * get(sampcount, munit, 0);
        if trted
            tcounts[φ] += 1 * get(sampcount, munit, 0);
        end
    end
end
