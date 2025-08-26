using TSCSMethods
using Test
using DataFrames
using Statistics

include("simulate_tscs.jl")

@testset "Coverage and Type I error" begin
    # Test parameters per specification
    F = 1:5    # Post-treatment periods
    L = -5:-1  # Pre-treatment periods
    n_seeds = 10  # Limited for PR testing (spec suggests 10 seeds)
    
    # Track coverage across seeds using null DGP (ATT = 0)
    coverage_results = []
    
    for seed in 1:n_seeds
        # Generate null data (DGP C with ATT = 0)
        df, treated_ids, t0_map = simulate_null(;
            N=600, T=200, treated_share=0.10, F=F, L=L,
            phi=0.0, rho=0.5, sigma_y=0.2, beta1=0.3, beta2=0.5, seed=seed
        )
        
        # Run TSCS pipeline
        timevary = Dict(:x1 => false, :x2 => true)
        model = makemodel(df, :t, :id, :gub, :y, [:x1, :x2], timevary, F, L)
        
        # Matching
        match!(model, df)
        
        # Estimate with higher iterations for coverage testing
        model.iterations = 400
        estimate!(model, df; dobayesfactor=false)
        
        # Check if 95% confidence intervals contain true effect (0)
        att_est = model.results.att
        ci_lower = model.results[!, Symbol("2.5%")]
        ci_upper = model.results[!, Symbol("97.5%")]
        
        # Coverage indicator for each F period (true effect = 0)
        coverage_per_f = (ci_lower .<= 0.0) .& (0.0 .<= ci_upper)
        push!(coverage_results, coverage_per_f)
    end
    
    # Compute empirical coverage across all seeds and F periods
    coverage_matrix = hcat(coverage_results...)  # Each column is a seed, rows are F periods
    empirical_coverage = mean(coverage_matrix, dims=2)[:, 1]  # Average across seeds for each F
    overall_coverage = mean(empirical_coverage)  # Overall coverage rate
    
    @info "Coverage Results" begin
        println("Empirical coverage per F period: ", round.(empirical_coverage, digits=3))
        println("Overall empirical coverage: ", round(overall_coverage, digits=3))
        println("Expected: ~0.95 for 95% CIs")
    end
    
    # PR-safe assertion per specification: coverage ∈ [0.90, 0.99]
    @test 0.90 <= overall_coverage <= 0.99 "95% CI empirical coverage should be in [0.90, 0.99]"
    
    # Additional check: most F periods should have reasonable coverage
    reasonable_coverage_count = sum(empirical_coverage .>= 0.85)
    @test reasonable_coverage_count >= length(F) ÷ 2 "Most F periods should have reasonable coverage"
end
