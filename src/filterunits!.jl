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
            m.matches[i].mus, m.matches[i].ranks, tt, tu, omap, m.F, m.ids,
            m.reference
        )
    end
end

function _filterunits!(mus, ranks, tt, tu, omap, F, ids, refval)
    rf = tt + refval
    if isnan(get(omap, (rf, tu), NaN)) | ismissing(get(omap, (rf, tu), missing))
        @views mus[:, :] .= false
        for ky in eachindex(ranks)
            ranks[ky] = Int[]
        end
        return 
    end

    for (j, f) in enumerate(F)

        ot = f + tt # outcome
        musj = @views mus[:, j]
        if isnan(get(omap, (ot, tu), NaN)) | ismissing(get(omap, (ot, tu), missing))
            musj .= false # whole column is trashed
            ranks[j] = Int[] # no ranks for trt obs at f
        else
            _filtermu!(musj, ranks[j], ot, refval, ids, omap)
        end
    end
end

function _filtermu!(musj, ranksj, ot, refval, ids, omap)
    for (l, mu) in enumerate(ids)
        if isnan(get(omap, (ot, mu), NaN)) | ismissing(get(omap, (ot, mu), missing)) | isnan(get(omap, (refval, mu), NaN)) | ismissing(get(omap, (refval, mu), missing))
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
