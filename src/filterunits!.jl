"""
filterunits!(m, omap)

Remove treated units and matches that do not have values defined in the outcome window. This could be modified to account for missingness easily.
This should probably be integrated into the match! procedure.
(keep separate for now)
"""
function filterunits!(m, omap)
    for (i, (tt, tu)) in enumerate(m.observations)
        for (j, f) in enumerate(m.F)
            ot = f + tt # outcome
            if isnan(get(omap, (ot, tu), NaN))
                m.matches[i].mus[:, j] .= false
            else
                for (l, mu) in enumerate(m.ids)
                    if isnan(get(omap, (ot, mu), NaN))
                        m.matches[i].mus[l, j] = false
                    end
                end
            end
        end
    end
end
