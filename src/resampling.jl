# Thread-safe buffer pool for sampling arrays
const SAMPLING_POOL = Ref{Union{Nothing, Channel{Vector{Int}}}}(nothing)

function __init_sampling_pool__()
    if SAMPLING_POOL[] === nothing
        pool = Channel{Vector{Int}}(Threads.nthreads() * 2)
        # Pre-populate with buffers
        for _ in 1:(Threads.nthreads() * 2)
            put!(pool, Vector{Int}(undef, 10000))  # Conservative max size
        end
        SAMPLING_POOL[] = pool
    end
end

function get_sampling_buffer(needed_size::Int)
    __init_sampling_pool__()
    
    if needed_size > 10000 || !isready(SAMPLING_POOL[])
        # Size too large or pool empty, allocate directly
        return Vector{Int}(undef, needed_size), false  # false = not pooled
    end
    
    buffer = take!(SAMPLING_POOL[])
    return resize!(buffer, needed_size), true  # true = pooled, needs return
end

function return_sampling_buffer(buffer::Vector{Int}, is_pooled::Bool)
    if is_pooled && length(buffer) <= 10000
        put!(SAMPLING_POOL[], buffer)
    end
end

"""
faster version of countmap, use when length is known
"""
function countmemb(itr, len::Int64)
  d = Dict{eltype(itr), Int64}()
  sizehint!(d, len)
  for val in itr
      d[val] = get(d, val, 0) + 1
  end
  return d
end

function sampleweights(observations, ids)
    sampcount = getsample(model)
    
    # Pre-allocate working arrays (could be pooled in future)
    matchweights = Vector{Int}(undef, length(ids))
    for (i, e) in enumerate(ids)
        matchweights[i] = get(sampcount, e, 0)
    end
    
    treatweights = Vector{Int}(undef, length(observations))
    for (i, (_, tunit)) in enumerate(observations)
        treatweights[i] = get(sampcount, tunit, 0)
    end
    
    return treatweights, matchweights
end

function treatedmap(observations)
    # Pre-size dictionary for better performance
    treatedmap = Dict{Int, Bool}()
    sizehint!(treatedmap, length(observations))
    for (_, tu) in observations
        treatedmap[tu] = true
    end
    return treatedmap
end

"""
not stratified
"""
function getsample(ids, treatedmap)
    # Use buffer pool to avoid allocation
    samp, is_pooled = get_sampling_buffer(length(ids))
    
    try
        # Fill with random samples
        sample!(ids, samp)
        executesample!(samp, treatedmap, ids)
        result = countmemb(samp, length(samp))
        return result
    finally
        # Always return buffer to pool
        return_sampling_buffer(samp, is_pooled)
    end
end

"""
not stratified
"""
function executesample!(samp, treatedmap, ids)

    while !checksample(samp, treatedmap)
        for j in eachindex(samp)
            samp[j] = sample(ids)
        end
    end

    return samp
end

"""
not stratified
"""
function checksample(samp, treatedmap)
    for unit in samp
        if get(treatedmap, unit, false)
            return true
        end
    end
    return false
end

# """
# stratified
# """
# function getsample(model::T) where T <: AbstractCICModelStratified
#     @unpack observations, strata, ids = model;
#     uniquestrata = sort(unique(strata));
#     stratdex = Dict(uniquestrata .=> 1:length(uniquestrata))
#     sbool = fill(false, length(uniquestrata))
    
#     treatedmap = Dict{Int, Tuple{Bool, Int}}();
#     for ((_, tunit), s) in zip(observations, strata)
#         treatedmap[tunit] = (true, s)
#     end

#     samp = sample(ids, length(ids));
#     executesample!(samp, sbool, treatedmap, ids, stratdex);
#     return countmemb(samp, length(samp))
# end

# """
# stratified
# """
# function executesample!(samp, sbool, treatedmap, ids, stratdex)

#     while !checksample(sbool, samp, treatedmap, stratdex)
#         for j in eachindex(samp)
#             samp[j] = sample(ids)
#         end
#     end

#     return samp
# end

# """
# stratified
# """
# function checksample(sbool, samp, treatedmap, stratdex)
#     for unit in samp
#         boolean, stratum = get(treatedmap, unit, (false, 1))
#         if boolean
#             sbool[get(stratdex, stratum, 0)] = true
#             if all(sbool)
#                 return true
#             end
#         end
#     end
#     return false
# end
