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
    model::AbstractCICModel, 
    dat::DataFrame;
    iterations::Union{Nothing, Int} = nothing,
    percentiles::Vector{Float64} = [0.025, 0.5, 0.975],
    overallestimate::Bool = false,
    dobayesfactor::Bool = true,
    dopvalue::Bool = false
)::Union{AbstractCICModel, Overall}

    # Input validation
    if nrow(dat) == 0
      throw(ArgumentError("Input data cannot be empty"))
    end
    
    if !isnothing(iterations) && iterations <= 0
      throw(ArgumentError("iterations must be positive"))
    end
    
    if !all(0 <= p <= 1 for p in percentiles)
      throw(ArgumentError("percentiles must be between 0 and 1"))
    end
    
    (; results, matches, observations, outcome, F, ids, reference, t, id) = model; 
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
        dat, Ys, Us,
        dobayesfactor, dopvalue
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
    else
        return model
    end
end

function _estimate!(
    results, matches, observations, outcome::Symbol,
    F, ids, reference, t, id, iterations, percentiles,
    dat, Ys, Us, dobayesfactor, dopvalue
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
    
    if dobayesfactor
        res[!, :bayesfactor] = fill(0.0, nrow(res))
        for (i, (tc, bot)) in enumerate(zip(Ys, eachrow(boots)))
            res[i, :bayesfactor] = bfactor(bot, tc)
        end
    end

    if dopvalue
        res[!, :pvalue] = fill(0.0, nrow(res))
        for (i, bot) in enumerate(eachrow(boots))
            res[i, :pvalue] = pvalue(bot)
        end
    end

    append!(results, res)

    # end lazy loop area for multiple outcomes
    return boots
end

"""
version for multiple outcomes

does not include p-values, bayesfactor
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
