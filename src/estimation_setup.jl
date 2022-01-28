# estimation_setup.jl

function processunits(model, dat)

    # for quick access to outcomes from unit & time
    outcomemap = getoutcomemap(dat, model);

    @unpack observations, matches, F, ids, reference = model;
    Fmin = minimum(F)
    obsnum = length(matches)

    # add total number of matches + treated unit
    # every relevant outcome (and reference)
    # over all treated units (and matches therein)
    for (k, mtch) in enumerate(matches)
        mtchsums[k] = sum(mtch.mus) + length(F)
    end

    # preallocate objects that hold the outcome-level
    # information, by treated obs
    # each element of these will be the set of relevant outcomes
    # for matches + treated unit for that treated observation
    Wos = Vector{Vector{Float64}}(undef, obsnum);
    Wrs = Vector{Vector{Float64}}(undef, obsnum);
    Tus = Vector{Vector{Int64}}(undef, obsnum);
    Mus = Vector{Vector{Int64}}(undef, obsnum);
    Fs = Vector{Vector{Int64}}(undef, obsnum);
    
    for (i, ms) in enumerate(mtchsums)
        Wos[i] = fill(0.0, ms)
        Wrs[i] = fill(0.0, ms)
        Tus[i] = fill(0, ms)
        Mus[i] = fill(0, ms)
        Fs[i] = fill(0, ms)
    end

    # populate these
    # essentially, use the match.mus information
    # to find the relevant outcomes and restructure them
    # for bootstrapping & estimation
    Threads.@threads for i in eachindex(observations)
        (tt, tu) = observations[i]
        mtch = matches[i]
        
        wo = Matrix{Union{Float64, Missing}}(missing, size(mtch.mus)) # outcomes
        wo2 = similar(wo) # reference

        wos = Wos[i]
        wrs = Wrs[i]
        tux = Tus[i]
        mux = Mus[i]
        fux = Fs[i]

        idx = findall(mtch.mus)
        matchnums = vec(sum(mtch.mus, dims = 1))

        unitstore!(
            wos, wrs, tux, mux, fux,
            wo, wo2,
            tt, tu,
            mtch.mus, ids, Fmin, outcomemap, matchnums, idx,
            reference
        )

    end
    return (Tus, Mus, Wos, Wrs, Fs)
end

# then do same for just treated, and append in
function unitstore!(
    wos, wrs, tux, mux, fux,
    wo, wo2,
    tt, tu,
    mus, ids, Fmin, outcomemap, matchnums, idx,
    reference
)
    for j in 1:size(mus, 2)
        ow = tt + j + Fmin - 1
        # treated unit
        # avals[j] += ((outcomemap[(tt-1, tu)] * -1.0) + (outcomemap[(ow, tu)] * 1.0))# * get(sampcount, tu, 0)
        for i in 1:size(mus, 1)
            if mus[i,j]
                id_i = @views ids[i];
                mn_j = @views matchnums[j];
                wo[i,j] = outcomemap[(ow, id_i)] * -inv(matchnums[j]) #* get(sampcount, ids[i], 0)# mult by wght and samp wght

                # matched unit
                wo2[i, j] = outcomemap[(tt+reference, id_i)] * inv(mn_j) #* get(sampcount, id_i, 0) # mult by wght and samp wght
            end
        end
    end

    # match units
    for (i, ix) in enumerate(idx)
        (r, c) = Tuple(ix)
        wos[i] = wo[ix]
        wrs[i] = wo2[ix]
        tux[i] = tu
        mux[i] = ids[r]
        fux[i] = c+Fmin-1
    end

    # treated unit
    trtx = length(idx)+1:length(idx)+length(F)
    for (g, tx) in enumerate(trtx)
        tux[tx] = tu
        mux[tx] = tu
        ow = tt + g + Fmin - 1
        fux[tx] = g + Fmin - 1
        wos[tx] = outcomemap[(ow, tu)] * 1.0 # treated unit weight, f in outcome window
        wrs[tx] = outcomemap[(tt+reference, tu)] * -1.0 # treated unit weight, reference
    end
   return wos, wrs, tux, mux, fux
end

function getoutcomemap(dat, model)
    outcomemap = Dict{Tuple{Int, Int}, Float64}();
    for r in eachrow(dat)
        outcomemap[(r[model.t], r[model.id])] = r[model.outcome]
    end

    return outcomemap
end

"""
        Fblock

Holds the relevant information for bootstrapping and estimation,
for a specific f in the outcome window (an f in a stratum when the model
is stratified).
"""
struct Fblock
    f::Int
    matchunits::Vector{Int}
    weightedoutcomes::Vector{Float64}
    weightedrefoutcomes::Vector{Float64}
    treatment::Vector{Bool}
end

"""
        makefblocks(subTus, subMus, subWos, subWrs, subFs)

Populate a set of fblocks to estimate and bootstrap the ATT.
"""
function makefblocks(subTus, subMus, subWos, subWrs, subFs)

    red = DataFrame(
        treatedunit = reduce(vcat, subTus),
        matchunit = reduce(vcat, subMus),
        woutcome = reduce(vcat, subWos),
        wref = reduce(vcat, subWrs),
        f = reduce(vcat, subFs)
    )

    red[!, :treatment] = red[!, :treatedunit] .== red[!, :matchunit]
    g1 = groupby(red, [:matchunit, :f]);
    g2 = @combine(
        g1,
        :woutcome = sum(:woutcome),
        :wref = sum(:wref),
        :treatment = any(:treatment)
    )

    G = groupby(g2, [:f], sort= true);

    fblocks = Vector{Fblock}(undef, length(G))

    for (i, gix) in enumerate(eachindex(G))
        g = G[gix]
        fblocks[i] = Fblock(
            gix[1],
            g.matchunit,
            g.woutcome,
            g.wref,
            g.treatment,
        )
    end
    return fblocks
end

"""
        stratifyinputs(X, s, strata)

Stratify the output from processunits().
"""
function stratifyinputs(X, s, strata)

    stratdex = strata .== s;

    return [@views(x[stratdex]) for x in X]
end