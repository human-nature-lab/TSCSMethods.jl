# sampling

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
    sampcount = getsample(model);
    matchweights = fill(0, length(ids));
    for (i, e) in enumerate(ids)
        matchweights[i] = get(sampcount, e, 0)
    end
    treatweights = fill(0, length(observations));
    for (i, (_, tunit)) in enumerate(observations)
        treatweights[i] = get(sampcount, tunit, 0)
    end
    return treatweights, matchweights
end

function treatedmap(observations)
    treatedmap = Dict{Int, Bool}();
    for (_, tu) in observations
        treatedmap[tu] = true
    end
    return treatedmap
end

"""
not stratified
"""
function getsample(ids, treatedmap)

    samp = sample(ids, length(ids));
    executesample!(samp, treatedmap, ids);
    return countmemb(samp, length(samp))
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

"""
stratified
"""
function getsample(model::T) where T <: AbstractCICModelStratified
    @unpack observations, strata, ids = model;
    uniquestrata = sort(unique(strata));
    stratdex = Dict(uniquestrata .=> 1:length(uniquestrata))
    sbool = fill(false, length(uniquestrata))
    
    treatedmap = Dict{Int, Tuple{Bool, Int}}();
    for ((_, tunit), s) in zip(observations, strata)
        treatedmap[tunit] = (true, s)
    end

    samp = sample(ids, length(ids));
    executesample!(samp, sbool, treatedmap, ids, stratdex);
    return countmemb(samp, length(samp))
end

"""
stratified
"""
function executesample!(samp, sbool, treatedmap, ids, stratdex)

    while !checksample(sbool, samp, treatedmap, stratdex)
        for j in eachindex(samp)
            samp[j] = sample(ids)
        end
    end

    return samp
end

"""
stratified
"""
function checksample(sbool, samp, treatedmap, stratdex)
    for unit in samp
        boolean, stratum = get(treatedmap, unit, (false, 1))
        if boolean
            sbool[get(stratdex, stratum, 0)] = true
            if all(sbool)
                return true
            end
        end
    end
    return false
end
