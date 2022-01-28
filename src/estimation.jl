"""
    bootinfo!(atts, boots; qtiles = [0.025, 0.5, 0.975])

Format the bootstrap matrix into the results dataframe. Assumes that att() has been executed.
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

"""
    estimate!(ccr::AbstractCICModel, dat; iterations = nothing)

Perform ATT estimation, with bootstrapped CIs.
"""
function estimate!(
    model::AbstractCICModel, dat;
    iterations = nothing, percentiles = [0.025, 0.5, 0.975]
)
    
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

Perform ATT estimation, with bootstrapped CIs.
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
        treatdex = treatedmap(observations);
        bootstrap!(
            multiboots[s], tcountmat, fblock_sub, ids, treatdex, iterations
        );

        multiatts[s] = fill(0.0, 31);
        tcounts = fill(0, 31);
        att!(multiatts[s], tcounts, fblock_sub)

        # add to dataframe
        for (r, e, f) in zip(eachrow(multiboots[s]), multiatts[s], F)
            push!(
                res,
                [
                    s, f,
                    e,
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

function setup_bootstrap(Flen, iterations)
    boots = zeros(Float64, Flen, iterations);
    tcountmat = zeros(Float64, Flen, iterations);
    return boots, tcountmat
end

function bootstrap!(boots, tcountmat, fblocks, ids, treatdex, iterations)
    @inbounds Threads.@threads for b in 1:iterations
        bootcol = @views boots[:, b]
        tcountcol = @views tcountmat[:, b]
        bootatt!(bootcol, tcountcol, fblocks, ids, treatdex);
    end
    return boots
end

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

## booting

# all fs, FOR A STRATUM NOT MULTIPLE
function bootatt!(atts, tcounts, fblocks, ids, treatdex)
    # b/c two methods w/ diff outputs
    # will need to separate anyway since atts will be different
    sampcount = getsample(ids, treatdex) # randomness
    _boot!(atts, tcounts, fblocks, sampcount)
end

function _boot!(atts, tcounts, fblocks, sampcount)
    for φ in 1:length(fblocks)
        @unpack matchunits, weightedoutcomes,
        weightedrefoutcomes, treatment = fblocks[φ]

        __boot!(
            atts, tcounts, φ,
            matchunits, weightedoutcomes,
            weightedrefoutcomes, treatment,
            sampcount
        )

        atts[φ] = atts[φ] * inv(tcounts[φ])
    end
end

function __boot!(
    atts, tcounts, φ,
    matchunits, weightedoutcomes,
    weightedrefoutcomes, treatments,
    sampcount
)
    for (munit, wo, wref, trted) in zip(
        matchunits, weightedoutcomes,
        weightedrefoutcomes, treatments
    )
        atts[φ] += (wo + wref) * get(sampcount, munit, 0);
        if trted
            tcounts[φ] += 1 * get(sampcount, munit, 0);
        end
    end
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
