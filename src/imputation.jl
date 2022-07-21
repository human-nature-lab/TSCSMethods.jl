"""
Impute Y(0).
"""
function impute_results(m, matches, dat, tvar; stratum = 1)

    res = m.results;
    oc = m.outcome;
    
    if m.stratifier != Symbol("")
        @subset!(res, :stratum .== 1)
    end

    odict = dict_mk(dat, m.outcome);

    sf = if m.stratifier != Symbol("")
        DataFrame(
            :timetreated => [e[1] for e in m.observations],
            :treatedunit => [e[2] for e in m.observations],
            :stratum => m.strata
        );
    else
        DataFrame(
            :timetreated => [e[1] for e in m.observations],
            :treatedunit => [e[2] for e in m.observations],
        );
    end
    
    mfo = leftjoin(matches, sf, on = [:timetreated, :treatedunit]);

    if m.stratifier != Symbol("")
        # toss out the excluded observations
        @subset!(mfo, :stratum .== stratum)
    end

    mfo = @transform(mfo, $(tvar) = :timetreated + :f);
    mfo[!, oc] .= 0.0;
    mfo[!, :match_values] .= 0.0;
    for r in eachrow(mfo)
        τ = r[:running]
        r[oc] = odict[(τ, r[:treatedunit])]
        mus = r[:matchunits]
        r[:match_values] = mean([odict[(τ, mu)] for mu in mus])
    end

    tobs = unique(mfo[!, [:timetreated, :treatedunit]]);
    mobs = unique(flatten(mfo, :matchunits)[!, [:timetreated, :matchunits]]);
    mobs[!, :τ] = mobs.timetreated .- 1;
    tobs[!, :τ] = tobs.timetreated .- 1;

    mobs[!, oc] .= 0.0;
    tobs[!, oc] .= 0.0;

    for u in eachrow(mobs)
        u[oc] = odict[(u[:τ], u[:matchunits])]
    end

    for u in eachrow(tobs)
        u[oc] = odict[(u[:τ], u[:treatedunit])]
    end

    mfosmy = @chain mfo begin
    groupby(:f)
    combine(
        oc => mean => oc,
        :match_values => mean => :match_values
    )
    end

    augres = leftjoin(res, mfosmy, on = :f);
    augres[!, :counter_post] = augres[!, oc] .- augres.att;

    # this needs to be fixed for abs
    mn, mx = sort([augres[!, oc] .- augres[!, Symbol("2.5%")], augres[!, oc] .- augres[!, Symbol("97.5%")]])

    augres[!, :counter_post_lwr] = mn
    augres[!, :counter_post_upr] = mx

    match_outcome_avg_pre = mean(mobs[!, m.outcome]);
    treated_outcome_avg_pre = mean(tobs[!, m.outcome]);

    augres[!, :pct] .= 0.0;
    augres[!, :pct_lo] .= 0.0;
    augres[!, :pct_hi] .= 0.0;

    for r in eachrow(augres)
        r[:pct] = change_pct(r[oc], r[:att])
        r[:pct_lo] = change_pct(r[oc], r[Symbol("2.5%")])
        r[:pct_hi] = change_pct(r[oc], r[Symbol("97.5%")])
    end

    # variability from day to day
    d2d = 100 .* diff(augres.match_values) .* inv.(augres.match_values[2:end]);
    augres[!, :pct_match_daytoday] = Vector{Union{Missing, Float64}}(undef, nrow(augres));
    augres[2:end, :pct_match_daytoday] = d2d;

    return augres, match_outcome_avg_pre, treated_outcome_avg_pre
end