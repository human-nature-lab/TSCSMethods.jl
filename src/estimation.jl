# estimation.jl

"""
    estimate!(
        ccr::AbstractCICModel, dat;
        iterations = nothing,
        percentiles = [0.025, 0.5, 0.975],
        overallestimate = false,
        bayesfactor = true
    )

Perform ATT estimation, with bootstrapped CIs.
"""
function estimate!(
    model::AbstractCICModel, dat;
    iterations = nothing,
    percentiles = [0.025, 0.5, 0.975],
    overallestimate = false,
    bayesfactor = true
)

    # import TSCSMethods:processunits,getoutcomemap,@unpack,unitstore!,setup_bootstrap,makefblocks,treatedmap,bootstrap!,att!,bootinfo!

    @unpack results, matches, observations, outcome, F, ids, reference, t, id = model; 
    modeliters = model.iterations;

    if !isnothing(iterations)
        @reset model.iterations = iterations;
    else
        iterations = modeliters
    end

    Ys, Us = unitcounts(model)

    boots = _estimate!(
        results, matches, observations, outcome,
        F, ids, reference, t, id, iterations, percentiles,
        dat, Ys, Us, bayesfactor
    );

    if overallestimate
        vb = vec(boots)

        return overall(
            att = mean(results.att),
            percentiles = quantile(vb, percentiles),
            bayesfactor = bfactor(vb, mean(results.treated)),
            ntreatedmean = mean(results.treated),
            pvalue = pvalue(vb)
        )
    end
end

import TSCSMethods:estimate!, _estimate!,bfactor,mean,pvalue,quantile,unitcounts,autobalance, make_groupindices, caliper, refine, meanbalance!, grandbalance!, checkbalances, checkwhile, refine

mutable struct Overall
    att::Float64
    percentiles::Vector{Float64}
    bayesfactor::Float64
    ntreatedmean::Float64
    pvalue::Float64
end

function overall(
    ; att = NaN, percentiles = Float64[], bayesfactor = NaN,
    ntreatedmean = 0.0, pvalue = NaN
)

    return Overall(
        att, percentiles, bayesfactor, ntreatedmean, pvalue
    )
end

function _estimate!(
    results, matches, observations, outcome::Symbol,
    F, ids, reference, t, id, iterations, percentiles,
    dat, Ys, Us, bayesfactor
)

    #=
    laziest way to do multiple outcomes:
        start a loop here over the vector of symbols
        turn X into Vector of typeof X, and append shit
        after first iter, make it add on to res
        add iter condition to nrow > 0 conditional
    =#

    X = processunits(
        matches, observations, outcome, F, ids, reference, t, id,
        dat
    );

    Flen = length(F);
    
    if nrow(results) > 0
        @reset results = DataFrame();
    end

    boots, tcountmat = setup_bootstrap(Flen, iterations)
    fblocks = makefblocks(X...)
    treatdex = treatedmap(observations);
    bootstrap!(boots, tcountmat, fblocks, ids, treatdex, iterations);

    atts = fill(0.0, Flen);
    tcounts = fill(0, Flen);
    att!(atts, tcounts, fblocks)

    res = DataFrame(f = F, att = atts)
    bootinfo!(res, boots; qtiles = percentiles)
    res[!, :treated] = Ys;
    res[!, :matches] = Us;
    
    if bayesfactor
        res[!, :bayesfactor] = fill(0.0, nrow(res))
        for (i, (tc, bot)) in enumerate(zip(Ys, eachrow(boots)))
            res[i, :bayesfactor] = bfactor(bot, tc)
        end
    end

    append!(results, res)

    # end lazy loop area for multiple outcomes
    return boots
end

"""
version for multiple outcomes
"""
function _estimate!(
    results, matches, observations, outcome::Vector{Symbol},
    F, ids, reference, t, id, iterations, percentiles,
    dat
)

    #=
    laziest way to do multiple outcomes:
        start a loop here over the vector of symbols
        turn X into Vector of typeof X, and append shit
        after first iter, make it add on to res
        add iter condition to nrow > 0 conditional
    =#

    res = DataFrame()

    for (ι, oc) in enumerate(outcome)

        X = processunits(
            matches, observations, oc, F, ids, reference, t, id,
            dat
        );
        
        if (nrow(results) > 0) & (ι == 1)
            @reset results = DataFrame();
        end

        boots, tcountmat = setup_bootstrap(length(F), iterations)
        fblocks = makefblocks(X...)
        treatdex = treatedmap(observations);
        bootstrap!(boots, tcountmat, fblocks, ids, treatdex, iterations);

        atts = fill(0.0, length(F));
        tcounts = fill(0, length(F));
        att!(atts, tcounts, fblocks)

        attsymb = Symbol("att" * "_" * string(oc))

        if ι == 1
            append!(
                res,
                DataFrame(
                    :f => F,
                    attsymb => atts
                )
            )
        else
            res[!, attsymb] = atts
        end

        bootinfo!(res, oc, boots; qtiles = percentiles)

    end
    
    append!(results, res)

    # end lazy loop area for multiple outcomes
end
