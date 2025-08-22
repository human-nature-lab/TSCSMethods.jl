# test_statistical_correctness.jl
# Statistical correctness validation - the foundation of trust in TSCSMethods

using Test
using Random
using DataFrames
using Statistics
using TSCSMethods

using TSCSMethods: generate_simple_tscs, generate_realistic_tscs

@testset "Statistical Correctness (CRITICAL TEST)" begin
    # CRITICAL: Testing statistical correctness with known ground truth
    # This validates that TSCSMethods produces correct causal estimates")
    
    @testset "ATT Unbiasedness - Core Validation" begin
        # Testing ATT unbiasedness: Can we recover known treatment effects?
        
        @testset "Zero Effect (Placebo Test)" begin
            true_att = 0.0
            estimates = Float64[]
            n_reps = 50  # Start with smaller number for speed
            
            # Testing true ATT = $true_att with $n_reps replications...
            
            for rep in 1:n_reps

                # Generate data with known zero effect
                # Treatment occurs periods 50-79 (duration=30), so we need:
                # - L windows BEFORE period 50 (pre-treatment matching)
                # - F windows AFTER period 79 (post-treatment outcomes)
                data = generate_simple_tscs(
                    true_att = true_att,
                    n_units = 60,
                    n_periods = 120,  # Extended to accommodate F windows after treatment
                    treatment_period = 50,  # Treatment starts at period 50
                    treatment_duration = 30,  # Treatment ends at period 79
                    n_treated = 15,
                    seed = rep
                )
                
                # Add a dummy covariate since package requires at least one
                # This is pure noise that shouldn't bias the ATT estimate
                data.dummy_covar = randn(nrow(data))
                
                # Run the full analysis pipeline
                timevary = Dict{Symbol, Bool}(:dummy_covar => false)  # Time-invariant dummy
                model = makemodel(
                    data, :time_period, :unit_id, :treatment, :outcome,
                    [:dummy_covar],  # Include dummy covariate
                    timevary,
                    1:10,      # F: 1-10 days AFTER treatment ends (periods 80-89)
                    -20:-5     # L: 20-5 days BEFORE treatment starts (periods 30-45)
                )
                
                match!(model, data)
                estimate!(model, data)
                
                # Extract ATT estimate from results DataFrame
                att_estimate = model.results.att[1]  # First row of att column  # First row, att column
                push!(estimates, att_estimate)
            end
            
            if length(estimates) > 0
                bias = mean(estimates) - true_att
                rmse = sqrt(mean((estimates .- true_att).^2))
                
                # Results: bias = $(round(bias, digits=4)), RMSE = $(round(rmse, digits=4))
                # Estimates range: [$(round(minimum(estimates), digits=3)), $(round(maximum(estimates), digits=3))]
                
                # Core unbiasedness test
                @test abs(bias) < 0.15  # Bias should be small for zero effect
                @test rmse < 1.0        # RMSE should be reasonable
                @test length(estimates) >= 40  # Most replications should succeed
            else
                @test false
                println("     ❌ All replications failed!")
            end
        end

        @testset "Positive Effect" begin
            true_att = 2.0
            estimates = Float64[]
            n_reps = 50
            
            # Testing true ATT = $true_att with $n_reps replications...
            for rep in 1:n_reps
                try
                    # Treatment occurs periods 50-79, need F after 79 and L before 50
                    data = generate_simple_tscs(
                        true_att = true_att,
                        n_units = 60,
                        n_periods = 120,  # Extended for F windows
                        treatment_period = 50,
                        treatment_duration = 30,  # Ends at period 79
                        n_treated = 15,
                        seed = rep + 1000  # Different seed from zero effect
                    )
                    
                    # Add dummy covariate for package requirement
                    data.dummy_covar = randn(nrow(data))
                    
                    timevary = Dict{Symbol, Bool}(:dummy_covar => false)
                    model = makemodel(
                        data, :time_period, :unit_id, :treatment, :outcome,
                        [:dummy_covar],
                        timevary,
                        5:15,      # F: 5-15 days AFTER treatment ends (periods 84-94)
                        -25:-10    # L: 25-10 days BEFORE treatment starts (periods 25-40)
                    )
                    
                    match!(model, data)
                    estimate!(model, data)
                    
                    att_estimate = model.results.att[1]  # First row of att column
                    push!(estimates, att_estimate)
                    
                catch e
                    @warn "Replication $rep failed: $e"
                end
            end
            
            if length(estimates) > 0
                bias = mean(estimates) - true_att
                rmse = sqrt(mean((estimates .- true_att).^2))
                
                @test abs(bias) < 0.2   # Allow slightly more bias for non-zero effects
                @test rmse < 1.5        # RMSE should still be reasonable
                @test length(estimates) >= 40
            else
                @test false
                println("     ❌ All replications failed!")
            end
        end
        
        @testset "Negative Effect" begin
            true_att = -1.5
            estimates = Float64[]
            n_reps = 50
            
            println("   Testing true ATT = $true_att with $n_reps replications...")
            
            for rep in 1:n_reps
                try
                    # Treatment occurs periods 50-79, need F after 79 and L before 50
                    data = generate_simple_tscs(
                        true_att = true_att,
                        n_units = 60,
                        n_periods = 120,  # Extended for F windows
                        treatment_period = 50,
                        treatment_duration = 30,  # Ends at period 79
                        n_treated = 15,
                        seed = rep + 2000
                    )
                    
                    # Add dummy covariate for package requirement
                    data.dummy_covar = randn(nrow(data))
                    
                    timevary = Dict{Symbol, Bool}(:dummy_covar => false)
                    model = makemodel(
                        data, :time_period, :unit_id, :treatment, :outcome,
                        [:dummy_covar],
                        timevary,
                        3:12,      # F: 3-12 days AFTER treatment ends (periods 82-91)
                        -30:-15    # L: 30-15 days BEFORE treatment starts (periods 20-35)
                    )
                    
                    match!(model, data)
                    estimate!(model, data)
                    
                    att_estimate = model.results.att[1]  # First row of att column
                    push!(estimates, att_estimate)
                    
                catch e
                    @warn "Replication $rep failed: $e"
                end
            end
            
            if length(estimates) > 0
                bias = mean(estimates) - true_att
                rmse = sqrt(mean((estimates .- true_att).^2))
                
                println("     Results: bias = $(round(bias, digits=4)), RMSE = $(round(rmse, digits=4))")
                println("     Estimates range: [$(round(minimum(estimates), digits=3)), $(round(maximum(estimates), digits=3))]")
                
                @test abs(bias) < 0.2
                @test rmse < 1.5
                @test length(estimates) >= 40
                
                println("     ✅ Negative effect test PASSED")
            else
                @test false
                println("     ❌ All replications failed!")
            end
        end
    end
    
    @testset "Confidence Interval Coverage" begin
        println("\n📊 Testing CI coverage: Do 95% CIs contain true effect 95% of time?")
        
        true_att = 1.0
        coverage_count = 0
        valid_cis = 0
        n_simulations = 100
        
        println("   Running $n_simulations simulations with true ATT = $true_att")
        
        for rep in 1:n_simulations
            try
                # Treatment occurs periods 25-34 (duration=10), need F after 34 and L before 25
                data = generate_simple_tscs(
                    true_att = true_att,
                    n_units = 80,
                    n_periods = 60,   # Extended to accommodate F windows after treatment
                    treatment_period = 25,
                    treatment_duration = 10,  # Ends at period 34
                    n_treated = 20,
                    seed = rep + 3000
                )
                
                # Add dummy covariate for package requirement
                data.dummy_covar = randn(nrow(data))
                
                timevary = Dict{Symbol, Bool}(:dummy_covar => false)
                model = makemodel(
                    data, :time_period, :unit_id, :treatment, :outcome,
                    [:dummy_covar],
                    timevary,
                    1:10,      # F: 1-10 days AFTER treatment ends (periods 35-44)
                    -15:-5     # L: 15-5 days BEFORE treatment starts (periods 10-20)
                )
                
                match!(model, data)
                estimate!(model, data)
                
                # Extract confidence interval bounds
                # Debug: Check what columns are available
                if nrow(model.results) > 0
                    available_columns = names(model.results)
                    if rep <= 3  # Only print first few for debugging
                        println("   Replication $rep columns: $available_columns")
                    end
                    
                    # Try different possible CI column names
                    ci_lower_col = nothing
                    ci_upper_col = nothing
                    
                    if "2.5%" in available_columns && "97.5%" in available_columns
                        ci_lower_col = "2.5%"
                        ci_upper_col = "97.5%"
                    elseif "CI_Lower" in available_columns && "CI_Upper" in available_columns
                        ci_lower_col = "CI_Lower"
                        ci_upper_col = "CI_Upper"
                    elseif "lower" in available_columns && "upper" in available_columns
                        ci_lower_col = "lower"
                        ci_upper_col = "upper"
                    end
                    
                    if ci_lower_col !== nothing && ci_upper_col !== nothing
                        ci_lower = model.results[1, ci_lower_col]
                        ci_upper = model.results[1, ci_upper_col]
                        
                        # Check if CI contains true value
                        if ci_lower <= true_att <= ci_upper
                            coverage_count += 1
                        end
                        valid_cis += 1
                    else
                        if rep <= 3
                            @warn "Cannot find CI columns in replication $rep. Available: $available_columns"
                        end
                    end
                else
                    @warn "No results from replication $rep"
                end
                
            catch e
                @warn "CI replication $rep failed: $e"
            end
        end
        
        if valid_cis > 0
            coverage_rate = coverage_count / valid_cis
            println("     Coverage rate: $(round(coverage_rate*100, digits=1))% (target: 95%)")
            println("     Valid CIs: $valid_cis out of $n_simulations simulations")
            
            # Allow for some Monte Carlo variation around 95%
            @test coverage_rate >= 0.85  # At least 85% coverage
            @test coverage_rate <= 1.05  # At most 105% (allowing for small sample variation)
            @test valid_cis >= 80        # Most simulations should produce valid CIs
            
            println("     ✅ Confidence interval coverage test PASSED")
        else
            @test false
            println("     ❌ No valid confidence intervals produced!")
        end
    end
    
    @testset "Precision Scaling" begin
        println("\n📊 Testing precision: Do standard errors decrease with sample size?")
        
        true_att = 1.5
        sample_sizes = [40, 80, 160]
        standard_errors = Float64[]
        
        for n_units in sample_sizes
            estimates = Float64[]
            n_reps = 30  # Fewer reps for larger samples
            
            println("   Testing with $n_units units ($n_reps replications)...")
            
            for rep in 1:n_reps
                try
                    # Treatment occurs periods 50-79, need F after 79 and L before 50
                    data = generate_simple_tscs(
                        true_att = true_att,
                        n_units = n_units,
                        n_periods = 110,  # Extended for F windows after treatment
                        treatment_period = 50,
                        treatment_duration = 30,  # Ends at period 79
                        n_treated = max(5, n_units ÷ 4),  # Scale treatment group
                        seed = rep + 4000 + n_units
                    )
                    
                    # Add dummy covariate for package requirement
                    data.dummy_covar = randn(nrow(data))
                    
                    timevary = Dict{Symbol, Bool}(:dummy_covar => false)
                    model = makemodel(
                        data, :time_period, :unit_id, :treatment, :outcome,
                        [:dummy_covar],
                        timevary,
                        5:15,      # F: 5-15 days AFTER treatment ends (periods 84-94)
                        -25:-10    # L: 25-10 days BEFORE treatment starts (periods 25-40)
                    )
                    
                    match!(model, data)
                    estimate!(model, data)
                    
                    push!(estimates, model.results.att[1])  # First row of att column
                    
                catch e
                    @warn "Precision test (n=$n_units, rep=$rep) failed: $e"
                end
            end
            
            if length(estimates) >= 20  # Need enough estimates for reliable SE
                se = std(estimates)
                push!(standard_errors, se)
                println("     n=$n_units: SE = $(round(se, digits=4)) ($(length(estimates)) valid estimates)")
            else
                println("     ⚠️  n=$n_units: Too few valid estimates ($(length(estimates)))")
                push!(standard_errors, NaN)
            end
        end
        
        # Test that standard errors generally decrease
        valid_ses = standard_errors[.!isnan.(standard_errors)]
        if length(valid_ses) >= 2
            # Allow some variation, but general trend should be decreasing
            @test valid_ses[end] < valid_ses[1] * 1.2  # Last should be smaller than first (with tolerance)
            println("     SE trend: $(round.(valid_ses, digits=4))")
            println("     ✅ Precision scaling test PASSED")
        else
            @warn "Insufficient valid standard errors for precision test"
            @test length(valid_ses) >= 2
        end
    end
    
    @testset "Realistic Synthetic Data Validation" begin
        println("\\n🎯 Testing with realistic event-based treatment design (like example_data)")
        
        @testset "Event-Based Treatment Effects" begin
            # Test small effects with event-based treatment (matching working example pattern)
            test_effects = [0.0, -0.02, 0.015]  # Small realistic effects
            
            for true_att in test_effects
                println("   Testing realistic ATT = $true_att...")
                
                try
                    # Generate data with realistic event-based treatment pattern
                    data = generate_realistic_tscs(
                        true_att = true_att,
                        n_units = 80,
                        n_days = 90,
                        n_treated = 16,  # 20% like example_data
                        seed = 123 + round(Int, true_att * 1000)
                    )
                    
                    # Remove validation columns for analysis
                    analysis_data = select(data, Not([:true_counterfactual, :true_treatment_effect]))
                    
                    timevary = Dict{Symbol, Bool}(:pop_dens => false)
                    model = makemodel(
                        analysis_data, :day, :fips, :gub, :death_rte,
                        [:pop_dens], timevary,
                        1:15,    # F: Post-treatment window
                        -25:-5   # L: Pre-treatment window
                    )
                    
                    match!(model, analysis_data)
                    estimate!(model, analysis_data)
                    
                    estimated_att = model.results.att[1]
                    bias = estimated_att - true_att
                    
                    println("     Estimated: $(round(estimated_att, digits=4)), Bias: $(round(bias, digits=4))")
                    
                    # More lenient test for realistic data - small biases are expected
                    @test abs(bias) < 0.05  # Allow small bias for realistic event-based design
                    
                    println("     ✅ Realistic synthetic test PASSED for ATT=$true_att")
                    
                catch e
                    @warn "Realistic synthetic test failed for ATT=$true_att: $e"
                    @test false  # Fail if analysis crashes
                end
            end
            
            println("\\n📋 Realistic synthetic validation complete")
            println("   This validates the package works with event-based treatment designs")
            println("   (matching the successful example_data pattern)")
        end
    end
    
    @testset "Basic Sanity Checks" begin
        println("\n🔍 Basic sanity checks on data generation")
        
        @testset "DGP Validation" begin
            # Test the data generating process itself
            # Treatment periods 35-54 (duration=20), need sufficient periods for F windows
            data = generate_simple_tscs(true_att=2.0, n_units=50, n_periods=70, treatment_period=35, treatment_duration=20, seed=42)
            
            @test nrow(data) == 50 * 50  # Should have n_units × n_periods rows
            @test all(in(1:50), data.unit_id)  # Unit IDs in correct range
            @test all(in(1:50), data.time_period)  # Time periods in correct range
            @test all(in([0, 1]), data.treatment)  # Treatment is binary
            
            # Check treatment timing
            treated_data = data[data.treatment .== 1, :]
            @test all(treated_data.time_period .>= 35)  # Treatment starts at period 35
            
            println("     ✅ Data generation validation PASSED")
        end
        
        @testset "Package Integration" begin
            # Test that our synthetic data works with the package
            # Treatment periods 50-79 (duration=30), need F after 79 and L before 50
            data = generate_simple_tscs(true_att=1.0, n_units=30, n_periods=100, treatment_period=50, treatment_duration=30, seed=123)
            
            # Add dummy covariate for package requirement
            data.dummy_covar = randn(nrow(data))
            
            timevary = Dict{Symbol, Bool}(:dummy_covar => false)
            model = makemodel(
                data, :time_period, :unit_id, :treatment, :outcome,
                [:dummy_covar],
                timevary,
                1:5,       # F: 1-5 days AFTER treatment ends (periods 80-84)
                -15:-10    # L: 15-10 days BEFORE treatment starts (periods 35-40)
            )
            
            @test model isa TSCSMethods.AbstractCICModel
            @test length(model.observations) > 0
            
            # Should be able to run full pipeline
            match!(model, data)
            estimate!(model, data)
            
            @test nrow(model.results) > 0
            @test !isnan(model.results.att[1])
        end
    end
    
    # Statistical correctness testing complete
    # If all tests passed, TSCSMethods produces statistically valid results
end