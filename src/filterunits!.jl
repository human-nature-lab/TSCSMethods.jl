"""
filterunits!(m, omap)

Remove treated units and matches that do not have values defined in the outcome window. This could be modified to account for missingness easily.
This should probably be integrated into the match! procedure.
(keep separate for now). (Updates ranks as well.)
"""
function filterunits!(m, omap)
    Threads.@threads for i in eachindex(m.observations)
        (tt, tu) = m.observations[i]

        _filterunits!(
            m.matches[i].mus, m.matches[i].ranks, tt, tu, omap, m.F, m.ids
        )
    end
end

function _filterunits!(mus, ranks, tt, tu, omap, F, ids)
    for (j, f) in enumerate(F)
        ot = f + tt # outcome
        musj = @views mus[:, j]
        if isnan(get(omap, (ot, tu), NaN))
            musj .= false # whole column is trashed
            ranks[j] = Int[] # no ranks for trt obs at f
        else
            _filtermu!(musj, ranks[j], ot, ids, omap)
        end
    end
end

function _filtermu!(musj, ranksj, ot, ids, omap)
    for (l, mu) in enumerate(ids)
        if isnan(get(omap, (ot, mu), NaN))
            musj[l] = false
            # ranks is in best-match order
            # find match to remove from ranks
            loc = findfirst(ranksj .== l)
            if !isnothing(loc)
                deleteat!(ranksj, loc)
            end
        end
    end
end
