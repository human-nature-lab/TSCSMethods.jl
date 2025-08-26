#!/usr/bin/env julia

using TSCSMethods
using DataFrames
include(joinpath(@__DIR__, "..", "test", "simulate_tscs.jl"))

function parse_range(s::AbstractString)
    # supports forms like "1:8" or "-30:-1"
    parts = split(s, ':')
    length(parts) == 2 || error("Invalid range: $s")
    a = parse(Int, parts[1]); b = parse(Int, parts[2])
    return a:b
end

function parse_args()
    args = Dict{String,String}()
    i = 1
    while i <= length(ARGS)
        if startswith(ARGS[i], "--")
            key = ARGS[i][3:end]
            i += 1
            i <= length(ARGS) || error("Missing value for --$key")
            args[key] = ARGS[i]
        else
            error("Unexpected arg: " * ARGS[i])
        end
        i += 1
    end
    return args
end

function to_json(x)
    # Minimal JSON encoder for simple types used here
    if x === nothing
        return "null"
    elseif x isa AbstractString
        return '"' * replace(x, '"' => "\\\"") * '"'
    elseif x isa Bool
        return string(x)
    elseif x isa Real
        return string(x)
    elseif x isa AbstractVector
        return "[" * join(map(to_json, x), ",") * "]"
    elseif x isa Dict
        parts = String[]
        for (k,v) in x
            push!(parts, '"' * String(k) * '"' * ":" * to_json(v))
        end
        return "{" * join(parts, ",") * "}"
    else
        return to_json(string(x))
    end
end

function main()
    a = parse_args()
    seeds = parse(Int, get(a, "seeds", "8"))
    N = parse(Int, get(a, "N", "1200"))
    T = parse(Int, get(a, "T", "220"))
    treated_share = parse(Float64, get(a, "treated-share", "0.12"))
    F = parse_range(get(a, "F", "1:5"))
    L = parse_range(get(a, "L", "-30:-1"))
    phi = parse(Float64, get(a, "phi", "0.0"))
    rho = parse(Float64, get(a, "rho", "0.5"))
    sigma_y = parse(Float64, get(a, "sigma-y", "0.12"))
    beta1 = parse(Float64, get(a, "beta1", "0.3"))
    beta2 = parse(Float64, get(a, "beta2", "0.5"))
    outpath = get(a, "out", "test/validation/seed_sweep.json")
    # Optional gates (not standard by default). Enable with --gate 1 or TSCS_GATES=1
    gate = get(a, "gate", get(ENV, "TSCS_GATES", "0")) == "1"
    mae_thresh = parse(Float64, get(a, "mae-thresh", "0.03"))
    mxe_thresh = parse(Float64, get(a, "mxe-thresh", "0.06"))

    delta = [0.00, 0.02, 0.04, 0.05, 0.05, 0.03, 0.02, 0.00][1:length(F)]

    results = Vector{Dict}(undef, seeds)
    maes = Float64[]
    mxes = Float64[]

    for s in 1:seeds
        df, _, _ = simulate_randomized(
            N=N, T=T, treated_share=treated_share, F=F, L=L, delta=delta,
            phi=phi, rho=rho, sigma_y=sigma_y, beta1=beta1, beta2=beta2, seed=100 + s
        )
        timevary = Dict(:x1 => false, :x2 => true)
        model = makemodel(df, :t, :id, :gub, :y, [:x1, :x2], timevary, F, L)
        match!(model, df)
        estimate!(model, df; dobayesfactor=false, iterations=100)
        est = collect(model.results.att)
        per_f = abs.(est .- delta)
        mae = mean(per_f)
        mxe = maximum(per_f)
        push!(maes, mae)
        push!(mxes, mxe)
        results[s] = Dict(
            "seed" => 100 + s,
            "mae" => mae,
            "mxe" => mxe,
            "per_f" => per_f
        )
    end

    summary = Dict(
        "mae_mean" => mean(maes),
        "mae_sd" => std(maes),
        "mxe_max" => maximum(mxes)
    )

    mkpath(dirname(outpath))
    open(outpath, "w") do io
        cfg = Dict(
            "seeds" => seeds, "N" => N, "T" => T, "treated_share" => treated_share,
            "F" => collect(F), "L" => collect(L), "phi" => phi, "rho" => rho,
            "sigma_y" => sigma_y, "beta1" => beta1, "beta2" => beta2
        )
        obj = Dict("config" => cfg, "results" => results, "summary" => summary)
        write(io, to_json(obj))
    end

    # By standard practice, report point-error metrics but do not gate by default
    if gate
        if summary["mae_mean"] > mae_thresh || summary["mxe_max"] > mxe_thresh
            @error "Seed sweep failed gates" summary
            exit(1)
        end
        println("Seed sweep OK (gated): ", summary)
    else
        println("Seed sweep report (no gates): ", summary)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
