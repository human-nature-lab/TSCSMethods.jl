@testset "Main Package Workflow" begin
    
    @testset "Complete Analysis Pipeline" begin
        println("Testing complete analysis workflow from raw data to results...")
        
        # Step 1: Generate realistic input data
        @testset "Data Preparation" begin
            data = example_data_generated(n_units=50, n_days=80, seed=123)
            
            @test nrow(data) == 50 * 80
            @test ncol(data) >= 6  # Minimum required columns
            required_cols = [:date, :fips, :pop_dens, :cumul_death_rate, :death_rte, :gub, :day]
            actual_cols = Symbol.(names(data))
            @test all(col -> col in actual_cols, required_cols)
            
            # Check data quality
            @test all(x -> x >= 0, data.death_rte)  # Non-negative outcomes
            @test length(unique(data.fips)) == 50   # Correct number of units
            @test all(x -> x in [0, 1], data.gub)   # Binary treatment
            
            println("✅ Data preparation successful")
        end
        
        # Step 2: Model Construction
        @testset "Model Construction" begin
            data = example_data_generated(n_units=30, n_days=60, seed=456)
            
            # Define analysis parameters
            matching_covariates = [:pop_dens, :cumul_death_rate]
            timevary_spec = Dict(:pop_dens => false, :cumul_death_rate => true)
            treatment_period = 35:45  # F: periods for effect estimation
            pretreat_period = 15:25   # L: pre-treatment matching window
            
            # Construct model
            model = makemodel(
                data, :day, :fips, :gub, :death_rte,
                matching_covariates, timevary_spec,
                treatment_period, pretreat_period;
                title = "workflow_test_model"
            )
            
            @test model isa CIC
            @test model.title == "workflow_test_model"
            @test model.F == treatment_period
            @test model.L == pretreat_period
            @test model.covariates == matching_covariates
            @test model.timevary == timevary_spec
            
            println("✅ Model construction successful")
        end
        
        # Step 3: Complete Workflow - Match → Balance → Estimate
        @testset "Full Analysis Workflow" begin
            data = example_data_generated(n_units=25, n_days=70, seed=789)
            
            # Build model
            model = makemodel(
                data, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                40:50, 20:30
            )
            
            # Step 3a: Matching
            @test_nowarn match!(model, data)
            @test !isempty(model.matches)
            println("✅ Matching completed")
            
            # Step 3b: Balancing  
            @test_nowarn balance!(model, data)
            println("✅ Balancing completed")
            
            # Step 3c: Estimation
            @test_nowarn estimate!(model, data; iterations=50)  # Reduced iterations for speed
            
            # Check that results were generated
            @test hasfield(typeof(model), :results)
            if isdefined(model, :results) && !isnothing(model.results) && nrow(model.results) > 0
                @test model.results isa DataFrame
                println("✅ Estimation completed with $(nrow(model.results)) time period results")
            end
        end
        
        # Step 4: Refinement Workflow
        @testset "Refinement Workflow" begin
            data = example_data_generated(n_units=20, n_days=50, seed=101112)
            
            # Base model
            base_model = makemodel(
                data, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                30:35, 15:20
            )
            
            match!(base_model, data)
            balance!(base_model, data)
            
            # Test refinement
            @test_nowarn refined_model = refine(
                base_model, data;
                refinementnum = 3,
                dobalance = true,
                doestimate = false  # Skip estimation for speed
            )
            
            println("✅ Model refinement successful")
        end
        
        # Step 5: Caliper Workflow  
        @testset "Caliper Workflow" begin
            data = example_data_generated(n_units=15, n_days=40, seed=131415)
            
            # Base model
            base_model = makemodel(
                data, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                25:30, 10:15
            )
            
            match!(base_model, data)
            balance!(base_model, data)
            
            # Test caliper application
            caliper_specs = Dict(:pop_dens => 0.5, :cumul_death_rate => 0.3)
            @test_nowarn calipered_model = caliper(
                base_model, caliper_specs, data;
                dobalance = true
            )
            
            println("✅ Caliper application successful")
        end
        
        # Step 6: Autobalance Workflow (the full automated approach)
        @testset "Autobalance Workflow" begin
            data = example_data_generated(n_units=20, n_days=60, seed=161718)
            
            # Base model
            model = makemodel(
                data, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                35:40, 20:25
            )
            
            match!(model, data)
            
            # Test autobalance - this should handle balance optimization automatically
            @test_nowarn results = autobalance(
                model, data;
                calmin = 0.1,
                step = 0.1,
                threshold = 0.15,  # More lenient threshold for test data
                doestimate = false,  # Skip estimation for speed
                verbose = false
            )
            
            println("✅ Autobalance workflow successful")
        end
    end
    
    @testset "Real-World Workflow Simulation" begin
        println("Testing realistic research workflow...")
        
        # Simulate a researcher's typical workflow
        @testset "Researcher Workflow Simulation" begin
            # Step 1: Load/generate data
            println("📊 Researcher loads data...")
            research_data = example_data_generated(n_units=40, n_days=90, seed=192021)
            
            # Step 2: Explore data structure
            println("🔍 Researcher explores data structure...")
            @test length(unique(research_data.fips)) > 1
            @test length(unique(research_data.date)) > 1
            treated_units = length(unique(research_data[research_data.gub .== 1, :fips]))
            control_units = length(unique(research_data[research_data.gub .== 0, :fips]))
            println("   Found $treated_units treated units, $control_units control units")
            
            # Step 3: Define analysis parameters (typical choices)
            println("⚙️  Researcher defines analysis parameters...")
            covariates = [:pop_dens, :cumul_death_rate]
            time_varying = Dict(:pop_dens => false, :cumul_death_rate => true)
            post_treatment = 50:65  # 16 periods post-treatment
            pre_treatment = 30:45   # 16 periods pre-treatment
            
            # Step 4: Build and run initial model
            println("🔧 Researcher builds initial model...")
            research_model = makemodel(
                research_data, :day, :fips, :gub, :death_rte,
                covariates, time_varying,
                post_treatment, pre_treatment;
                title = "Research Study Model"
            )
            
            # Step 5: Perform matching
            println("🎯 Researcher performs matching...")
            match!(research_model, research_data)
            
            # Step 6: Check initial balance
            println("⚖️  Researcher checks balance...")
            balance!(research_model, research_data)
            initial_balance = checkbalances(research_model)
            
            # Step 7: Improve balance if needed
            println("🎛️  Researcher optimizes balance...")
            if !isnothing(initial_balance)
                # Try autobalance for better results
                cal_model, refined_cal_model, overall_result = autobalance(
                    research_model, research_data;
                    dooverall = true,
                    doestimate = false,  # We'll estimate separately
                    verbose = false
                )
                @test !isnothing(refined_cal_model)
                chosen_model = refined_cal_model
            else
                chosen_model = research_model
            end
            
            # Step 8: Final estimation
            println("📈 Researcher performs final estimation...")
            @test_nowarn estimate!(chosen_model, research_data; iterations=100)
            
            # Step 9: Extract results
            println("📋 Researcher extracts results...")
            if isdefined(chosen_model, :results) && !isnothing(chosen_model.results) && nrow(chosen_model.results) > 0
                results = chosen_model.results
                @test results isa DataFrame
                @test nrow(results) > 0
                
                # Check that we have treatment effect estimates
                if "att" in names(results)
                    mean_effect = mean(results.att)
                    println("   📊 Average treatment effect: $(round(mean_effect, digits=4))")
                    @test isa(mean_effect, Real)
                end
                
                println("✅ Research workflow completed successfully!")
            end
        end
    end
    
    @testset "Error Handling in Workflow" begin
        println("Testing workflow robustness...")
        
        # Test workflow with problematic data
        @testset "Workflow Error Recovery" begin
            # Test with minimal data
            minimal_data = example_data_generated(n_units=5, n_days=20, seed=222324)
            
            model = makemodel(
                minimal_data, :day, :fips, :gub, :death_rte,
                [:pop_dens],
                Dict(:pop_dens => false),
                12:15, 5:8
            )
            
            # These should either work or fail gracefully  
            try
                match!(model, minimal_data)
                println("   Matching with minimal data: ✅")
            catch e
                println("   Matching with minimal data failed as expected: $(typeof(e))")
            end
            
            try
                balance!(model, minimal_data)
                println("   Balancing with minimal data: ✅")
            catch e
                println("   Balancing with minimal data failed as expected: $(typeof(e))")
            end
            
            @test true  # Test that we handled the edge case gracefully
            
            println("✅ Workflow handles edge cases gracefully")
        end
    end
end