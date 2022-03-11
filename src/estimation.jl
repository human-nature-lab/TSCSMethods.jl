# estimation.jl

"""
        att!(atts, tcounts, fblocks)

Calculate the att for each f, for a the set of treated units and
matches contained in the fblocks.
"""
function att!(atts, tcounts, fblocks)
    for φ in 1:length(fblocks)
        @unpack matchunits, weightedoutcomes,
        weightedrefoutcomes, treatment = fblocks[φ]

        __att!(
            atts, tcounts, φ,
            weightedoutcomes,
            weightedrefoutcomes, treatment
        )

        atts[φ] = atts[φ] * inv(tcounts[φ])
    end
end

function __att!(
    atts, tcounts, φ,
    weightedoutcomes,
    weightedrefoutcomes, treatments,
)
    for (wo, wref, trted) in zip(
        weightedoutcomes,
        weightedrefoutcomes, treatments
    )
        atts[φ] += (wo + wref);
        if trted
            tcounts[φ] += 1;
        end
    end
end

"""
    estimate!(ccr::AbstractCICModel, dat; iterations = nothing)

Perform ATT estimation, with bootstrapped CIs.
"""
function estimate!(
    model::AbstractCICModel, dat;
    iterations = nothing,
    percentiles = [0.025, 0.5, 0.975]
)

    import TSCSMethods:processunits,getoutcomemap,@unpack,unitstore!,setup_bootstrap,makefblocks,treatedmap,bootstrap!,att!,bootinfo!,applyunitcounts!

    X = processunits(model, dat);
    @unpack observations, matches, F, ids = model;
    Flen = length(F);

    if !isnothing(iterations)
        @reset model.iterations = iterations;
    else
        iterations = model.iterations
    end
    
    if nrow(model.results) > 0
        @reset model.results = DataFrame();
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

    append!(model.results, res)

    applyunitcounts!(model)

    return model
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

# utilities

"""
    bootinfo!(res, boots; qtiles = [0.025, 0.5, 0.975])

Format the bootstrap matrix into the results dataframe. Assumes that att()
has already been added to res.
"""
function bootinfo!(res, boots; qtiles = [0.025, 0.5, 0.975])
  qnmes = Vector{Symbol}();
  for q in qtiles
    res[!, :mean] = Vector{Float64}(undef, nrow(res))
    qn = Symbol(string(q * 100) * "%");
    push!(qnmes, qn)
    res[!, qn] = Vector{Float64}(undef, nrow(res))
  end
  
  for (c, r) in enumerate(eachrow(res))
    r[qnmes] = quantile(boots[c, :], qtiles)
    r[:mean] = mean(boots[c, :])
  end
  return res
end

function applyunitcounts!(model)
  
  Ys, Us = unitcounts(model)

  res = model.results;
  res[!, :treated] = zeros(Int, nrow(res))
  res[!, :matches] = zeros(Int, nrow(res))

  strat = any(
    [
        typeof(model) == x for x in [
            CICStratified, RefinedCICStratified, RefinedCaliperCICStratified
        ]
    ]
  );

  if !strat
    Yd = Dict(collect(1:length(model.F)) .=> Ys)
    Ud = Dict(collect(1:length(model.F)) .=> Us)

    for (j, f) in enumerate(res.f)
      φ = f - minimum(model.F) + 1
      res[j, :treated] = Ys[φ]
      res[j, :matches] = Us[φ]
    end
  else
    for s in res.stratum
      Yd = Dict(collect(1:length(model.F)) .=> Ys[s])
      Ud = Dict(collect(1:length(model.F)) .=> Us[s])
      
      for (j, f, s) in zip(eachindex(res.f), res.f, res.stratum)
        φ = f - minimum(model.F) + 1
        res[j, :treated] = Ys[s][φ]
        res[j, :matches] = Us[s][φ]
      end
    end
  end
  return model
end

function matchprocess(mfo, dat; ovars = [deathoutcome, caseoutcome])

    vn = VariableNames();

    mfo[!, :matchnum] .= 0
    for (i, e) in enumerate(mfo.matchunits)
        mfo.matchnum[i] = length(e)
    end

    mfol = flatten(mfo, [:matchunits, :ranks]);

    tfo = unique(mfol[!, [:timetreated, :treatedunit, :f]]);
    tfo[!, :matchunits] = tfo[!, :treatedunit]
    tfo[!, :ranks] = fill(0, nrow(tfo))
    tfo[!, :matchnum] = fill(0, nrow(tfo))

    mfol = vcat(mfol, tfo);
    sort!(mfol, [:timetreated, :treatedunit, :f, :matchunits]);

    
    mfol = unique(mfol[!, [:timetreated, :treatedunit, :matchunits, :ranks, :matchnum]])
    rename!(mfol, :matchunits => :matchunit, :ranks => :rank)
    
    mfol[!, :treated] = mfol[!, :matchunit] .== mfol[!, :treatedunit]
    
    # only works with whole range present
    for o in ovars
        mfol[!, o] = Vector{Vector{Float64}}(undef, nrow(mfol))
    end

    for r in eachrow(mfol)
        mu = r[:matchunit]
        tt = r[:timetreated]
    
        c1 = (dat[!, vn.id] .== mu);
        ct = (dat[!, vn.t] .>= tt - 30) .& (dat[!, vn.t] .<=  tt + 80);
    
        for o in ovars
            r[o] = dat[c1 .& ct, o]
        end
    end
    
    return mfol
end

"""
quick_att(matchseries;
    outcomes = [:death_rte, :case_rte, :deaths, :cases],
    F = 10:40,
    ttt = false,
    tm1 = 30
)

Quickly calculate the att for specified outcomes based on the output of the 
matchprocess(). Does not calculate confidence intervals.

Calculates both the individual-level unit effects, and the average effects.

"""
function quick_att(
    matchseries;
    outcomes = [:death_rte, :case_rte, :deaths, :cases],
    F = 10:40,
    ttt = false,
    tm1 = 30
)

    # (280:350)[31]

    attdat = select(
        matchseries,
        :timetreated, :treatedunit, :matchunit;
    );

    gg = unique(attdat[!, [:treatedunit, :timetreated]]);
    
    for ot in outcomes
        attdat[!, ot] = Vector{Vector{Float64}}(undef, nrow(attdat))
        gg[!, ot] = Vector{Vector{Float64}}(undef, nrow(gg))
    end

    # r = matchseries[1, :];
    for (i, r) in enumerate(eachrow(matchseries))
        for ot in outcomes
            # change within units: Y_(t+F) - Y_(t-1)
            attdat[i, ot] = [r[ot][φ + 31] - r[ot][tm1] for φ in F]
        end
    end

    # assign negative values to the matched units
    # for the difference in changes
    asub = @views attdat[attdat.treatedunit .!= attdat.matchunit, :];
    asub[!, outcomes] = asub[!, outcomes] .* -1;

    # take the average of the matches to a single treated unit
    # (this will expode the time series vars
    # should be order preserving)

    # unit difference in changes?
    # matches first x calcualtion
    # attdat.death_rte[3][1] + mean([attdat.death_rte[i][1] for i in [1,2,4,5,6]])

    # avg. over match units
    x = @chain attdat begin
        @transform(:treated = :treatedunit .== :matchunit)
        groupby([:timetreated, :treatedunit, :treated])
        combine(
            [ot => mean => ot for ot in outcomes]...
        )
    end
    
    
    ln = Int(nrow(x) * inv(length(F)))

    # add the φ values, in order
    x[!, :φ] = reduce(vcat, fill(collect(F), ln))

    
    # SUM NOT AVERAGE....
    # individual unit effect
    x = @chain x begin
        groupby([:timetreated, :treatedunit, :φ])
        combine(
            [ot => sum => ot for ot in outcomes]...
        )
    end

    func = if ttt
        sum
    else mean
    end

    # atts
    oa = @chain x begin
        groupby(:φ)
        combine(
            [ot => func => ot for ot in outcomes]...
        )
    end

    return oa, x
end