# estimation_setup.jl

"""
        getoutcomemap(dat, model)

Get dictionary, from specified unit and time to outcome.
"""
function getoutcomemap(dat, model)

    outcomemap = Dict{Tuple{Int, Int}, Float64}()

    # there shouldn't be any missing outcomes actually pulled
    # outcomemap = if Missing <: eltype(dat[!, model.outcome])
    #     Dict{Tuple{Int, Int}, Union{Float64, Missing}}()
    # else
    #     Dict{Tuple{Int, Int}, Float64}()
    # end
    for r in eachrow(dat)
        if !ismissing(r[model.outcome])
            outcomemap[(r[model.t], r[model.id])] = r[model.outcome]
        end
    end

    return outcomemap
end

function processunits(model, dat)

    # for quick access to outcomes from unit & time
    outcomemap = getoutcomemap(dat, model);

    @unpack observations, matches, F, ids, reference = model;
    Fmin = minimum(F)
    obsnum = length(matches)

    # add total number of matches + treated unit
    # every relevant outcome (and reference)
    # over all treated units (and matches therein)
    mtchsums = Vector{Int}(undef, length(matches));
    for (k, mtch) in enumerate(matches)

        # number of columns with any matches
        # gives the number of times the treated unit exists
        Φ = 0
        for col in eachcol(mtch.mus); if any(col); Φ += 1 end end

        mtchsums[k] = sum(mtch.mus) + Φ
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
    
    Threads.@threads for i in eachindex(mtchsums)
        ms = mtchsums[i]
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
        
        #wo = Matrix{Union{Float64, Missing}}(missing, size(mtch.mus)) # outcomes
        #wo2 = similar(wo) # reference

        wos = Wos[i]
        wrs = Wrs[i]
        tux = Tus[i]
        mux = Mus[i]
        fux = Fs[i]

        matchnums = vec(sum(mtch.mus, dims = 1))

        unitstore!(
            wos, wrs, tux, mux, fux,
            tt, tu,
            mtch.mus, ids, Fmin, outcomemap, matchnums,
            reference
        )
    end
    return (Tus, Mus, Wos, Wrs, Fs)
end

"""

Use the mu matrix to generate vectors of weighted outcomes
and unit information, as a preprocessing step for estimation.
"""
function unitstore!(
    wos, wrs, tux, mux, fux,
    tt, tu,
    mus, ids, Fmin, outcomemap, matchnums,
    reference
)

    # k is sort of a linear index for the matrix
    # counts true elements in mu (matches to treated), and
    # an (additional) an element for every column (for the treated unit)
    k = 0
    for (j, col) in enumerate(eachcol(mus))
        if any(col)
            ow = tt + j + Fmin - 1
            k += 1
            # treated unit
            wos[k] = outcomemap[(ow, tu)] * 1.0 # treated unit weight, f in outcome window
            wrs[k] = outcomemap[(tt+reference, tu)] * -1.0 # treated unit weight, reference
            tux[k] = tu
            mux[k] = tu
            fux[k] = j + Fmin - 1
            for (i, e) in enumerate(col)
                if e
                    k += 1
                    # matches to treated
                    wos[k] = outcomemap[(ow, ids[i])] * -inv(matchnums[j])
                    wrs[k] = outcomemap[(tt+reference, ids[i])] * inv(matchnums[j])
                    tux[k] = tu
                    mux[k] = ids[i]
                    fux[k] = j + Fmin - 1
                end
            end
        end
    end

   return wos, wrs, tux, mux, fux
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
