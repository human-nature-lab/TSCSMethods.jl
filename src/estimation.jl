# estimation.jl

"""
    estimate!(ccr::AbstractCICModel, dat; iterations = nothing)

Perform ATT estimation, with bootstrapped CIs.
"""
function estimate!(
    model::AbstractCICModel, dat;
    iterations = nothing,
    percentiles = [0.025, 0.5, 0.975]
)

    # import TSCSMethods:processunits,getoutcomemap,@unpack,unitstore!,setup_bootstrap,makefblocks,treatedmap,bootstrap!,att!,bootinfo!,applyunitcounts!

    @unpack results, matches, observations, outcome, F, ids, reference, t, id = model; 
    modeliters = model.iters;

    if !isnothing(iterations)
        @reset model.iterations = iterations;
    else
        iterations = modeliters
    end


    _estimate!(
        results, matches, observations, outcome::Symbol,
        F, ids, reference, t, id, iterations, percentiles,
        dat
    );

    applyunitcounts!(model)

    return model
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
            res = DataFrame(
                :f => F,
                attsymb => atts
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

    multiboots = Dict{Int, Matrix{Float64}}();
    multiatts = Dict{Int, Vector{Float64}}();
    
    res = DataFrame(
        stratum = Int[], f = Int[], att = Float64[], mean = Float64[]
    )

    for q in percentiles
        qn = Symbol(string(q * 100) * "%");
        res[!, qn] = Float64[]
    end

    X = processunits(model, dat);
    @unpack observations, matches, F, ids = model;
    Flen = length(F);

    if !isnothing(iterations)
        @reset model.iterations = iterations;
    else
        iterations = model.iterations
    end
    
    if (nrow(model.results) > 0) | length(names(model.results)) > 0
        @reset model.results = DataFrame();
    end

    for s in sort(unique(model.strata))
        Xsub = stratifyinputs(X, s, model.strata)
        multiboots[s], tcountmat = setup_bootstrap(Flen, iterations)
        fblock_sub = makefblocks(Xsub...)
        obsub = @views observations[model.strata .== s]
        treatdex = treatedmap(obsub);
        bootstrap!(
            multiboots[s], tcountmat, fblock_sub, ids, treatdex, iterations
        );

        multiatts[s] = fill(0.0, length(F));
        tcounts = fill(0, length(F));
        att!(multiatts[s], tcounts, fblock_sub)

        # add to dataframe
        for (r, e, f) in zip(eachrow(multiboots[s]), multiatts[s], F)
            push!(
                res,
                [
                    s, f, e,
                    mean(r),
                    (quantile(r, percentiles))...
                ]
            )
        end
    end

    append!(model.results, res)
    
    applyunitcounts!(model)

    return model
end
