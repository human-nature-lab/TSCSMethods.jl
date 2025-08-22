"""
    estimate!(ccr::AbstractCICModelStratified, dat; iterations = nothing)

Perform ATT stratified estimation, with bootstrapped CIs.
"""
function estimate!(
    model::AbstractCICModelStratified, dat;
    iterations = nothing,
    percentiles = [0.025, 0.5, 0.975],
    overallestimate = false,
    dobayesfactor = false,
    dopvalue = false

)

    (; results, matches, observations, strata, outcome, F, ids, reference, t, id) = model; 
    modeliters = model.iterations;

    if !isnothing(iterations)
        @reset model.iterations = iterations;
    else
        iterations = modeliters
    end

    Ys, Us = unitcounts(model)

    multiboots = Dict{Int, Matrix{Float64}}();
    multiatts = Dict{Int, Vector{Float64}}();

    _estimate_strat!(
        multiatts, multiboots,
        results, matches, observations, strata, outcome,
        F, ids, reference, t, id, iterations, percentiles,
        dat, Ys, Us, dobayesfactor, dopvalue
    )

    stra = sort(unique(strata))
    
    if overallestimate
        oe = overall(
            att = fill(NaN, length(stra)),
            percentiles = Vector{Vector{Float64}}(undef, length(stra)),
            bayesfactor = fill(NaN, length(stra)),
            ntreatedmean = fill(NaN, length(stra)),
            pvalue = fill(NaN, length(stra)),
            stratum = stra
        )

        for (l, s) in enumerate(stra)

            mbs = vec(multiboots[s])

            oe.att[l] = mean(multiatts[s])
            oe.percentiles[l] = quantile(mbs, percentiles)
            # CHANGE
            oe.bayesfactor[l] = bfactor(mbs, mean(results.treated[results.stratum .== s]))
            oe.ntreatedmean[l] = mean(results.treated[results.stratum .== s])
            oe.pvalue[l] = pvalue(mbs)

        end
        return oe
    end
end

function _estimate_strat!(
    multiatts, multiboots,
    results, matches, observations, strata, outcome::Vector{Symbol},
    F, ids, reference, t, id, iterations, percentiles,
    dat, Ys, Us,
    dobayesfactor, dopvalue
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

                if dobayesfactor
                    res_add[!, :bayesfactor] = [bfactor(r, tc)]
                end

                if dopvalue
                    res_add[!, :pvalue] = [pvalue(r)]
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
    dat, Ys, Us,
    dobayesfactor, dopvalue
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

            if dobayesfactor
                res_add[!, :bayesfactor] = [bfactor(r, tc)]
            end

            if dopvalue
                res_add[!, :pvalue] = [pvalue(r)]
            end

            append!(
                res,
                res_add
            )
        end
    end

    append!(results, res)
end
