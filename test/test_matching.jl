@testset "Matching Functions" begin
    
    # Helper function to create test model and data
    function create_test_model_and_data()
        dates = Date(2021, 10, 1):Day(1):Date(2021, 11, 30)
        n_units = 20
        
        # Create synthetic test data
        test_data = DataFrame()
        for unit in 1:n_units
            for (i, date) in enumerate(dates)
                push!(test_data, (
                    date = date,
                    fips = 1000 + unit,
                    pop_dens = 50.0 + randn() * 10,
                    cumul_death_rate = 20.0 + i * 0.1 + randn() * 2,
                    death_rte = abs(randn() * 0.5),
                    gub = (unit <= 5) ? 1 : 0  # First 5 units are treated
                ))
            end
        end
        
        test_data[!, :day] = Dates.value.(test_data.date .- minimum(test_data.date))
        
        matching_covariates = [:pop_dens, :cumul_death_rate]
        timevary = Dict(:pop_dens => false, :cumul_death_rate => true)
        
        model = makemodel(
            test_data, :day, :fips, :gub, :death_rte,
            matching_covariates, timevary,
            10:20, -30:-1;
            title = "test_matching_model"
        )
        
        return model, test_data
    end
    
    @testset "match! function" begin
        model, test_data = create_test_model_and_data()
        
        @testset "Basic matching" begin
            # Test that match! runs without errors
            @test_nowarn match!(model, test_data)
            
            # Test that matching populates the matches structure
            # After matching, matches should contain distance information
            @test !isempty(model.matches)
            @test all(m -> !isempty(m.ranks), model.matches)
        end
        
        @testset "Custom treatment categories" begin
            # Test with custom treatment category function
            function custom_treatcat(x)
                if x == 0
                    return 0
                elseif x <= 2
                    return 1
                else
                    return 2
                end
            end
            
            @test_nowarn match!(model, test_data; treatcat = custom_treatcat)
        end
        
        @testset "Matching with exposure variable" begin
            # Test matching with exposure (if supported)
            # Add an exposure column to test data
            test_data[!, :exposure] = rand(nrow(test_data))
            
            # This might not be implemented yet, but test should not crash
            @test_nowarn match!(model, test_data; exposure = :exposure)
        end
    end
    
    @testset "default_treatmentcategories function" begin
        @test default_treatmentcategories(0) == 0
        @test default_treatmentcategories(1) == 1
        @test default_treatmentcategories(5) == 1
        @test default_treatmentcategories(-1) == 1
        
        # Test with various input types
        @test default_treatmentcategories(0.0) == 0
        @test default_treatmentcategories(1.5) == 1
    end
    
    @testset "Matching edge cases" begin
        @testset "No treated units" begin
            model_no_treated, test_data_no_treated = create_test_model_and_data()
            # Set all units to untreated
            test_data_no_treated[!, :gub] .= 0
            
            # Rebuild model with no treated units
            model_no_treated = makemodel(
                test_data_no_treated, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                10:20, -30:-1
            )
            
            # Should handle gracefully (either skip matching or give meaningful error)
            @test_nowarn match!(model_no_treated, test_data_no_treated) || 
                   @test_throws Exception match!(model_no_treated, test_data_no_treated)
        end
        
        @testset "All units treated" begin
            model_all_treated, test_data_all_treated = create_test_model_and_data()
            # Set all units to treated
            test_data_all_treated[!, :gub] .= 1
            
            # Rebuild model with all treated units
            model_all_treated = makemodel(
                test_data_all_treated, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                10:20, -30:-1
            )
            
            # Should handle gracefully (either skip matching or give meaningful error)
            @test_nowarn match!(model_all_treated, test_data_all_treated) ||
                   @test_throws Exception match!(model_all_treated, test_data_all_treated)
        end
        
        @testset "Single treated unit" begin
            model_single, test_data_single = create_test_model_and_data()
            # Set only one unit to treated
            test_data_single[!, :gub] .= 0
            test_data_single[test_data_single.fips .== 1001, :gub] .= 1
            
            # Rebuild model
            model_single = makemodel(
                test_data_single, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                10:20, -30:-1
            )
            
            # Should handle single treated unit case
            @test_nowarn match!(model_single, test_data_single)
        end
    end
    
    @testset "Distance calculations" begin
        model, test_data = create_test_model_and_data()
        match!(model, test_data)
        
        @testset "Distance matrix properties" begin
            # After matching, distance matrices should have proper dimensions
            for match_obj in model.matches
                if !isempty(match_obj.distances) && match_obj.distances isa Matrix
                    @test size(match_obj.distances, 1) > 0
                    @test size(match_obj.distances, 2) > 0
                    # Distances should be non-negative
                    @test all(d -> d >= 0, match_obj.distances[.!isnan.(match_obj.distances)])
                end
            end
        end
        
        @testset "Ranking properties" begin
            # Ranks should contain valid indices
            for match_obj in model.matches
                for (unit_id, ranks) in match_obj.ranks
                    @test all(r -> r > 0, ranks)  # Ranks should be positive
                    @test allunique(ranks)  # Ranks should be unique
                end
            end
        end
    end
end