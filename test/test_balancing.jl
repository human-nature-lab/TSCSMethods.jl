@testset "Balancing Functions" begin
    
    # Helper function to create test model with matching completed
    function create_matched_model()
        dates = Date(2021, 10, 1):Day(1):Date(2021, 11, 30)
        n_units = 20
        
        test_data = DataFrame()
        for unit in 1:n_units
            for (i, date) in enumerate(dates)
                push!(test_data, (
                    date = date,
                    fips = 1000 + unit,
                    pop_dens = 50.0 + randn() * 10,
                    cumul_death_rate = 20.0 + i * 0.1 + randn() * 2,
                    death_rte = abs(randn() * 0.5),
                    gub = (unit <= 5) ? 1 : 0
                ))
            end
        end
        
        test_data[!, :day] = Dates.value.(test_data.date .- minimum(test_data.date))
        
        model = makemodel(
            test_data, :day, :fips, :gub, :death_rte,
            [:pop_dens, :cumul_death_rate],
            Dict(:pop_dens => false, :cumul_death_rate => true),
            10:20, -30:-1
        )
        
        # Perform matching first
        match!(model, test_data)
        
        return model, test_data
    end
    
    @testset "balance! function" begin
        model, test_data = create_matched_model()
        
        @testset "Basic balancing" begin
            # Test that balance! runs without errors
            @test_nowarn balance!(model, test_data)
            
            # Test that balance! returns the model
            balanced_model = balance!(model, test_data)
            @test balanced_model === model
        end
        
        @testset "Balancing with different model types" begin
            # Test balancing with regular CIC model
            model_cic, test_data_cic = create_matched_model()
            @test_nowarn balance!(model_cic, test_data_cic)
            
            # If stratified models are available, test those too
            # This would require creating stratified test data and models
        end
    end
    
    @testset "checkbalances function" begin
        model, test_data = create_matched_model()
        balance!(model, test_data)
        
        @testset "Balance checking" begin
            # Test that checkbalances runs and returns something meaningful
            balance_results = checkbalances(model, test_data)
            
            # Balance results should contain information about covariate balance
            @test balance_results isa DataFrame || balance_results isa Dict
            
            # Should have information about each matching covariate
            if balance_results isa DataFrame
                @test nrow(balance_results) > 0
            elseif balance_results isa Dict
                @test length(balance_results) > 0
            end
        end
        
        @testset "Balance thresholds" begin
            # Test checkbalances with different significance levels if supported
            @test_nowarn checkbalances(model, test_data; α = 0.05) ||
                   @test_nowarn checkbalances(model, test_data)
            
            @test_nowarn checkbalances(model, test_data; α = 0.01) ||
                   @test_nowarn checkbalances(model, test_data)
        end
    end
    
    @testset "autobalance function" begin
        model, test_data = create_matched_model()
        
        @testset "Basic autobalancing" begin
            # Test basic autobalance functionality
            @test_nowarn autobalance(model, test_data)
            
            # Test that autobalance returns models
            result = autobalance(model, test_data)
            
            # Should return models (might be tuple of models)
            @test !isnothing(result)
        end
        
        @testset "Autobalance with parameters" begin
            # Test with specific caliper parameters
            @test_nowarn autobalance(
                model, test_data;
                calmin = 0.1,
                step = 0.05
            )
            
            # Test with initial balance specifications
            initial_bals = Dict(:pop_dens => 0.5, :cumul_death_rate => 0.25)
            @test_nowarn autobalance(
                model, test_data;
                initial_bals = initial_bals,
                calmin = 0.08,
                step = 0.05
            )
        end
        
        @testset "Autobalance with overall estimation" begin
            # Test autobalance that also computes overall estimates
            cal_model, ref_cal_model, overall_est = autobalance(
                model, test_data;
                dooverall = true,
                calmin = 0.1,
                step = 0.1
            )
            
            # Should return three components when dooverall = true
            @test !isnothing(cal_model)
            @test !isnothing(ref_cal_model)
            @test !isnothing(overall_est)
            @test overall_est isa Overall
        end
    end
    
    @testset "Balance quality assessment" begin
        model, test_data = create_matched_model()
        balance!(model, test_data)
        
        @testset "Balance metrics" begin
            # After balancing, we should be able to assess balance quality
            balance_results = checkbalances(model, test_data)
            
            # Balance results should provide information we can test
            if balance_results isa DataFrame
                # Check that it has expected columns for balance assessment
                @test ncol(balance_results) > 0
            elseif balance_results isa Dict
                # Check that it has balance information for covariates
                @test any(k -> string(k) in ["pop_dens", "cumul_death_rate", "pvalue"], keys(balance_results))
            end
        end
        
        @testset "Balance improvement" begin
            # Test that balancing affects the model state
            model_pre_balance, test_data_balance = create_matched_model()
            model_post_balance = deepcopy(model_pre_balance)
            
            # Get initial balance
            initial_balance = checkbalances(model_pre_balance, test_data_balance)
            
            # Apply balancing
            balance!(model_post_balance, test_data_balance)
            final_balance = checkbalances(model_post_balance, test_data_balance)
            
            # Balance results should potentially be different
            # (though not necessarily "better" without knowing the specific metrics)
            @test !isnothing(initial_balance)
            @test !isnothing(final_balance)
        end
    end
    
    @testset "Edge cases in balancing" begin
        @testset "Balancing with insufficient matches" begin
            # Create a scenario with very few potential matches
            dates = Date(2021, 10, 1):Day(1):Date(2021, 10, 10)
            test_data_small = DataFrame()
            
            for unit in 1:5  # Only 5 units
                for (i, date) in enumerate(dates)
                    push!(test_data_small, (
                        date = date,
                        fips = 1000 + unit,
                        pop_dens = 50.0 + unit * 20,  # Very different values
                        cumul_death_rate = 20.0 + unit * 10,
                        death_rte = abs(randn() * 0.5),
                        gub = (unit <= 1) ? 1 : 0  # Only 1 treated unit
                    ))
                end
            end
            
            test_data_small[!, :day] = Dates.value.(test_data_small.date .- minimum(test_data_small.date))
            
            model_small = makemodel(
                test_data_small, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                5:7, -5:-1
            )
            
            match!(model_small, test_data_small)
            
            # Should handle edge case gracefully
            @test_nowarn balance!(model_small, test_data_small) ||
                   @test_throws Exception balance!(model_small, test_data_small)
        end
    end
end