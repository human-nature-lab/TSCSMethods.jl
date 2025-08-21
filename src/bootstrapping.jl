# boostrapping.jl

function setup_bootstrap(Flen, iterations)
    boots = zeros(Float64, Flen, iterations);
    tcountmat = zeros(Float64, Flen, iterations);
    return boots, tcountmat
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

        atts[φ] = atts[φ] * inv(tcounts[φ])
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
