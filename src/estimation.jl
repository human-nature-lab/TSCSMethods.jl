# estimation.jl

"""
    estimate!(ccr::AbstractCICModel, dat; iterations = nothing)

Perform ATT estimation, with bootstrapped CIs.
"""
function estimate!(
    model::AbstractCICModel, dat;
    iterations = nothing,
    percentiles = [0.025, 0.5, 0.975],
    bootout = false,
)

    # import TSCSMethods:processunits,getoutcomemap,@unpack,unitstore!,setup_bootstrap,makefblocks,treatedmap,bootstrap!,att!,bootinfo!,applyunitcounts!

    @unpack results, matches, observations, outcome, F, ids, reference, t, id = model; 
    modeliters = model.iterations;

    if !isnothing(iterations)
        @reset model.iterations = iterations;
    else
        iterations = modeliters
    end


    boots = _estimate!(
        results, matches, observations, outcome,
        F, ids, reference, t, id, iterations, percentiles,
        dat
    );

    applyunitcounts!(model)

    if bootout
        return boots
    end
end

function _estimate!(
    results, matches, observations, outcome::Symbol,
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

"""
    estimate!(ccr::AbstractCICModelStratified, dat; iterations = nothing)

Perform ATT stratified estimation, with bootstrapped CIs.
"""
function estimate!(
    model::AbstractCICModelStratified, dat;
    iterations = nothing, percentiles = [0.025, 0.5, 0.975]
)

    @unpack results, matches, observations, strata, outcome, F, ids, reference, t, id = model; 
    modeliters = model.iterations;

    if !isnothing(iterations)
        @reset model.iterations = iterations;
    else
        iterations = modeliters
    end

    multiboots = Dict{Int, Matrix{Float64}}();
    multiatts = Dict{Int, Vector{Float64}}();

    _estimate_strat!(
        multiatts, multiboots,
        results, matches, observations, strata, outcome,
        F, ids, reference, t, id, iterations, percentiles,
        dat
    )
    
    applyunitcounts!(model)

    return model
end

function _estimate_strat!(
    multiatts, multiboots,
    results, matches, observations, strata, outcome::Vector{Symbol},
    F, ids, reference, t, id, iterations, percentiles,
    dat
)
    
    if (nrow(results) > 0) | length(names(results)) > 0
        @reset results = DataFrame();
    end

    reses = DataFrame()

    for (ι, oc) in enumerate(outcome)

        X = processunits(
            matches, observations, oc, F, ids, reference, t, id,
            dat
        );
            
        res = DataFrame()
        Flen = length(F);

        for s in sort(unique(strata))
            Xsub = stratifyinputs(X, s, strata)
            multiboots[s], tcountmat = setup_bootstrap(Flen, iterations)
            fblock_sub = makefblocks(Xsub...)
            obsub = @views observations[strata .== s]
            treatdex = treatedmap(obsub);
            bootstrap!(
                multiboots[s], tcountmat, fblock_sub, ids, treatdex, iterations
            );

            multiatts[s] = fill(0.0, length(F));
            tcounts = fill(0, length(F));
            att!(multiatts[s], tcounts, fblock_sub)

            attsymb = Symbol("att" * "_" * string(oc))
            prefix = string(oc) * "_"
            barname = Symbol(prefix * "mean")
            losymb = Symbol(prefix * string(percentiles[1] * 100) * "%");
            medsymb = Symbol(prefix * string(percentiles[2] * 100) * "%");
            hisymb = Symbol(prefix * string(percentiles[3] * 100) * "%");

            # add to dataframe
            for (r, e, f) in zip(eachrow(multiboots[s]), multiatts[s], F)
                (lov, miv, hiv) = quantile(r, percentiles)
                append!(
                    res,
                    DataFrame(
                        :stratum => s,
                        :f => f,
                        attsymb => e,
                        barname => mean(r),
                        losymb => lov,
                        medsymb => miv,
                        hisymb => hiv
                    )
                )
            end
        end

        if ι == 1
            append!(reses, res)
        else
            reses = leftjoin(reses, res, on = [:stratum, :f])
        end

    end
    
    append!(results, reses)
end

function _estimate_strat!(
    multiatts, multiboots,
    results, matches, observations, strata, outcome::Symbol,
    F, ids, reference, t, id, iterations, percentiles,
    dat
)
    
    if (nrow(results) > 0) | length(names(results)) > 0
        @reset results = DataFrame();
    end


    X = processunits(
        matches, observations, outcome, F, ids, reference, t, id,
        dat
    );
        
    res = DataFrame()
    Flen = length(F);

    for s in sort(unique(strata))
        Xsub = stratifyinputs(X, s, strata)
        multiboots[s], tcountmat = setup_bootstrap(Flen, iterations)
        fblock_sub = makefblocks(Xsub...)
        obsub = @views observations[strata .== s]
        treatdex = treatedmap(obsub);
        bootstrap!(
            multiboots[s], tcountmat, fblock_sub, ids, treatdex, iterations
        );

        multiatts[s] = fill(0.0, length(F));
        tcounts = fill(0, length(F));
        att!(multiatts[s], tcounts, fblock_sub)

        attsymb = :att
        barname = :mean
        losymb = Symbol(string(percentiles[1] * 100) * "%");
        medsymb = Symbol(string(percentiles[2] * 100) * "%");
        hisymb = Symbol(string(percentiles[3] * 100) * "%");

        # add to dataframe
        for (r, e, f) in zip(eachrow(multiboots[s]), multiatts[s], F)
            (lov, miv, hiv) = quantile(r, percentiles)
            append!(
                res,
                DataFrame(
                    :stratum => s,
                    :f => f,
                    attsymb => e,
                    barname => mean(r),
                    losymb => lov,
                    medsymb => miv,
                    hisymb => hiv
                )
            )
        end
    end

    append!(results, res)
end