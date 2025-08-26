#!/usr/bin/env julia

using TSCSMethods
using DataFrames
using Statistics

include(joinpath(@__DIR__, "..", "test", "simulate_tscs.jl"))

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

function run_bias_analysis(seeds::Int, iterations::Int)
    F = 1:5    # Post-treatment periods  
    L = -5:-1  # Pre-treatment periods
    delta = [0.01, 0.02, 0.05, 0.03, 0.01]  # Known true effects
    
    bias_results = []
    
    for seed in 1:seeds
        # Generate randomized data (DGP A) with known delta
        df, _, _ = simulate_randomized(;
            N=600, T=200, treated_share=0.10, F=F, L=L, delta=delta,
            phi=0.0, rho=0.5, sigma_y=0.2, beta1=0.3, beta2=0.5, seed=seed
        )
        
        # Run TSCS pipeline
        timevary = Dict(:x1 => false, :x2 => true)
        model = makemodel(df, :t, :id, :gub, :y, [:x1, :x2], timevary, F, L)
        
        match!(model, df)
        estimate!(model, df; dobayesfactor=false, iterations=iterations)
        
        est = collect(model.results.att)
        bias = est .- delta
        rmse_per_f = sqrt.(bias.^2)
        
        push!(bias_results, Dict(
            "seed" => seed,
            "bias" => bias,
            "rmse" => rmse_per_f,
            "mae" => mean(abs.(bias)),
            "overall_rmse" => sqrt(mean(bias.^2))
        ))
    end
    
    # Compute summary statistics
    all_bias = hcat([r["bias"] for r in bias_results]...)
    mean_bias = mean(all_bias, dims=2)[:, 1]
    rmse_bias = sqrt.(mean(all_bias.^2, dims=2))[:, 1]
    
    return Dict(
        "results" => bias_results,
        "summary" => Dict(
            "mean_bias_per_f" => mean_bias,
            "rmse_per_f" => rmse_bias,
            "overall_mae" => mean([r["mae"] for r in bias_results]),
            "overall_rmse" => mean([r["overall_rmse"] for r in bias_results])
        )
    )
end

function run_coverage_analysis(seeds::Int, iterations::Int)
    F = 1:5    # Post-treatment periods
    L = -5:-1  # Pre-treatment periods
    
    coverage_results = []
    
    for seed in 1:seeds
        # Generate null data (DGP C) with ATT = 0
        df, _, _ = simulate_null(;
            N=600, T=200, treated_share=0.10, F=F, L=L,
            phi=0.0, rho=0.5, sigma_y=0.2, beta1=0.3, beta2=0.5, seed=seed
        )
        
        # Run TSCS pipeline
        timevary = Dict(:x1 => false, :x2 => true)
        model = makemodel(df, :t, :id, :gub, :y, [:x1, :x2], timevary, F, L)
        
        match!(model, df)
        estimate!(model, df; dobayesfactor=false, iterations=iterations)
        
        # Check 95% CI coverage of true effect (0)
        ci_lower = model.results[!, Symbol("2.5%")]
        ci_upper = model.results[!, Symbol("97.5%")]
        coverage_per_f = (ci_lower .<= 0.0) .& (0.0 .<= ci_upper)
        
        push!(coverage_results, Dict(
            "seed" => seed,
            "coverage" => coverage_per_f,
            "overall_coverage" => mean(coverage_per_f)
        ))
    end
    
    # Compute empirical coverage
    all_coverage = hcat([r["coverage"] for r in coverage_results]...)
    empirical_coverage_per_f = mean(all_coverage, dims=2)[:, 1]
    overall_coverage = mean(empirical_coverage_per_f)
    
    return Dict(
        "results" => coverage_results,
        "summary" => Dict(
            "coverage_per_f" => empirical_coverage_per_f,
            "overall_coverage" => overall_coverage,
            "target" => 0.95
        )
    )
end

function main()
    a = parse_args()
    seeds = parse(Int, get(a, "seeds", "20"))
    iterations = parse(Int, get(a, "iterations", "400"))
    outpath = get(a, "out", "test/validation/bias_coverage.json")
    
    println("Running bias analysis with $seeds seeds and $iterations iterations...")
    bias_section = run_bias_analysis(seeds, iterations)
    
    println("Running coverage analysis with $seeds seeds and $iterations iterations...")
    coverage_section = run_coverage_analysis(seeds, iterations)
    
    # Combine results
    results = Dict(
        "bias" => bias_section,
        "coverage" => coverage_section,
        "config" => Dict(
            "seeds" => seeds,
            "iterations" => iterations
        )
    )
    
    # Write output
    mkpath(dirname(outpath))
    open(outpath, "w") do io
        write(io, to_json(results))
    end
    
    # Report summary (standard practice: gate on coverage only)
    coverage_rate = coverage_section["summary"]["overall_coverage"]
    println("Coverage analysis - Overall coverage: ", round(coverage_rate, digits=4))

    # Coverage gates near nominal 95%
    if coverage_rate < 0.93 || coverage_rate > 0.97
        @error "Coverage outside acceptable [0.93, 0.97]" coverage_rate
        exit(1)
    end

    println("Bias and coverage analysis completed successfully (coverage within [0.93, 0.97])")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
