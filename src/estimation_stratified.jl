"""
    estimate!(ccr::AbstractCICModelStratified, dat; iterations = nothing)

Perform ATT stratified estimation, with bootstrapped CIs.
"""
function estimate!(
    model::AbstractCICModelStratified, dat;
    iterations = nothing,
    percentiles = [0.025, 0.5, 0.975],
    overall = false,
    bayesfactor = true
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
        dat, bayesfactor
    )

    overalls = Dict{Int, Tuple{Float64, Vector{Float64}}}()
    overalls_bf = Dict{Int, Float64}()
    if overall
        for s in sort(unique(strata))
            overalls_bf[s] = bfactor(vec(multiboots[s]), mean(results.treated))
            overalls[s] = (
                mean(multiatts[s]),
                quantile(vec(multiboots[s]), percentiles)
            )
        end
        return overalls, overalls_bf
    end
end

function _estimate_strat!(
    multiatts, multiboots,
    results, matches, observations, strata, outcome::Vector{Symbol},
    F, ids, reference, t, id, iterations, percentiles,
    dat, bayesfactor
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

        Ys, Us = unitcounts(model)

        for s in sort(unique(strata))
            Ysub = Ys[s]
            Usub = Us[s]
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
            for (r, e, f, tc, mc) in zip(
                eachrow(multiboots[s]), multiatts[s], F,
                Ysub, Usub
            )
                (lov, miv, hiv) = quantile(r, percentiles)

                res_add = DataFrame(
                    :stratum => s,
                    :f => f,
                    attsymb => e,
                    barname => mean(r),
                    losymb => lov,
                    medsymb => miv,
                    hisymb => hiv,
                    :treated => tc,
                    :matches => mc
                )

                if bayesfactor
                    res_add[!, :bayesfactor] = [bfactor(r, tc)]
                end

                append!(
                    res,
                    res_add
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
    dat, bayesfactor
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

    Ys, Us = unitcounts(model)

    for s in sort(unique(strata))
        Ysub = Ys[s]
        Usub = Us[s]
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
        for (r, e, f, tc, mc) in zip(
            eachrow(multiboots[s]), multiatts[s], F, Ysub, Usub
        )
            (lov, miv, hiv) = quantile(r, percentiles)

            res_add = DataFrame(
                :stratum => s,
                :f => f,
                attsymb => e,
                barname => mean(r),
                losymb => lov,
                medsymb => miv,
                hisymb => hiv,
                :treated => tc,
                :matches => mc
            )

            if bayesfactor
                res_add[!, :bayesfactor] = [bfactor(r, tc)]
            end

            append!(
                res,
                res_add
            )
        end
    end

    append!(results, res)
end
