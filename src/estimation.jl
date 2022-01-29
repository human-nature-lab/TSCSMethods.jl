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
        treatdex = treatedmap(observations);
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
