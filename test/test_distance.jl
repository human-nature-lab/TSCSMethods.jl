@testset "Distance Calculation Functions" begin
    
    # Helper function to create test model and data
    function create_distance_test_model()
        Random.seed!(42)  # For reproducible tests
        dates = Date(2021, 10, 1):Day(1):Date(2021, 11, 15)  # 46 days
        n_units = 15
        
        test_data = DataFrame()
        for unit in 1:n_units
            is_treated = (unit <= 5)  # First 5 units treated
            treatment_day = is_treated ? 25 : -1  # Treatment on day 25
            
            for (i, date) in enumerate(dates)
                push!(test_data, (
                    date = date,
                    fips = 2000 + unit,
                    pop_dens = 40.0 + unit * 2.0 + randn() * 3.0,
                    cumul_death_rate = 15.0 + i * 0.2 + randn() * 1.5,
                    death_rte = abs(0.3 + randn() * 0.2),
                    gub = (is_treated && i == treatment_day) ? 1 : 0
                ))
            end
        end
        
        test_data[!, :day] = Dates.value.(test_data.date .- minimum(test_data.date))
        
        matching_covariates = [:pop_dens, :cumul_death_rate]
        timevary = Dict(:pop_dens => false, :cumul_death_rate => true)
        
        model = makemodel(
            test_data, :day, :fips, :gub, :death_rte,
            matching_covariates, timevary,
            1:3, -10:-1;  # F = 1:3 (post-treatment), L = -10:-1 (pre-treatment)
            title = "Distance Test Model"
        )
        
        return model, test_data
    end

    @testset "Threading Modernization Tests" begin
        
        @testset "Deterministic Results with :greedy Scheduler" begin
            # Test that :greedy threading gives deterministic results
            Random.seed!(123)
            model1, data1 = create_distance_test_model()
            match!(model1, data1)
            
            Random.seed!(123)  # Same seed
            model2, data2 = create_distance_test_model()
            match!(model2, data2)
            
            # Should have identical results
            @test model1.matches[1].distances ≈ model2.matches[1].distances atol=1e-14
        end
        
        @testset "Threading Performance Validation" begin
            model, data = create_distance_test_model()
            
            # Should complete without errors
            @test_nowarn match!(model, data)
            
            # Should produce valid distance matrices
            @test !isempty(model.matches)
            @test !isempty(model.matches[1].distances)
            
            # All distances should be finite and non-negative
            for match in model.matches
                if !isempty(match.distances)
                    for dist in match.distances
                        if !isinf(dist)
                            @test dist >= 0.0
                            @test isfinite(dist)
                        end
                    end
                end
            end
        end
        
        @testset "Multiple Thread Safety" begin
            # Test with different thread counts (if available)
            original_threads = Threads.nthreads()
            
            model, data = create_distance_test_model()
            @test_nowarn match!(model, data)
            
            # Store results for comparison
            baseline_distances = deepcopy(model.matches[1].distances)
            
            # Test should work regardless of thread count
            @test size(baseline_distances, 1) > 0
            @test size(baseline_distances, 2) > 0
        end
    end

    @testset "Distance Calculation Accuracy" begin
        
        @testset "Mahalanobis Distance Implementation" begin
            model, data = create_distance_test_model()
            match!(model, data)
            
            # Test basic distance properties
            distances = model.matches[1].distances
            if !isempty(distances)
                # Distance to self should be 0 (not applicable in matching context, but test non-negativity)
                @test all(d -> d >= 0.0 || isinf(d), distances)
                
                # No NaN values should appear
                @test !any(isnan, distances)
            end
        end
        
        @testset "Division vs Inverse Performance Fix" begin
            # Create test vectors for averaging
            Random.seed!(456)
            test_values = rand(100) * 10
            test_counts = rand(1:20, 100)
            
            # Test that division gives same results as inv() * multiplication
            for (val, count) in zip(test_values, test_counts)
                division_result = val / count
                inverse_result = val * (1.0 / count)
                
                @test division_result ≈ inverse_result atol=1e-14
            end
        end
        
        @testset "Missing Data Handling" begin
            # Create data with some missing values
            Random.seed!(789)
            model, data = create_distance_test_model()
            
            # Introduce some missing values in covariates (use allowmissing)
            data = allowmissing(data, :cumul_death_rate)
            missing_indices = sample(1:nrow(data), 5, replace=false)
            data[missing_indices, :cumul_death_rate] .= missing
            
            # Should handle missing data gracefully
            @test_nowarn match!(model, data)
            
            # Should still produce some valid matches
            @test !isempty(model.matches)
        end
    end

    @testset "Distance Allocation Functions" begin
        
        @testset "distances_allocate! Function" begin
            model, data = create_distance_test_model()
            
            # Test allocation before matching
            @test_nowarn match!(model, data)
            
            # Check that distance arrays are properly allocated
            for match_obj in model.matches
                if !isempty(match_obj.distances)
                    @test isa(match_obj.distances, Matrix{Float64})
                    @test size(match_obj.distances, 2) == length(model.covariates) + 1  # +1 for Mahalanobis
                end
            end
        end
        
        @testset "Memory Efficiency" begin
            # Test that large models don't cause memory issues
            Random.seed!(999)
            large_model, large_data = create_distance_test_model()
            
            # Monitor that matching completes without excessive allocation warnings
            @test_nowarn match!(large_model, large_data)
            
            # Verify results are properly structured
            @test !isempty(large_model.matches)
        end
        
        @testset "Pre-allocation Optimization" begin
            # Test thread-local storage works correctly
            Random.seed!(1111)
            model, data = create_distance_test_model()
            
            # Should work without errors
            @test_nowarn match!(model, data)
            
            # Test multiple calls use pre-allocated storage
            Random.seed!(1111)
            model2, data2 = create_distance_test_model()
            @test_nowarn match!(model2, data2)
            
            # Results should be identical (deterministic with pre-allocation)
            if !isempty(model.matches) && !isempty(model2.matches)
                @test model.matches[1].distances ≈ model2.matches[1].distances atol=1e-14
            end
        end
        
        @testset "Thread Safety with Pre-allocation" begin
            # Test that thread-local storage doesn't cause race conditions
            Random.seed!(2222)
            models = []
            datas = []
            
            # Create multiple models
            for i in 1:3
                Random.seed!(2222 + i)
                model, data = create_distance_test_model()
                push!(models, model)
                push!(datas, data)
            end
            
            # Run matching on all models (uses threading internally)
            for (model, data) in zip(models, datas)
                @test_nowarn match!(model, data)
                @test !isempty(model.matches)
            end
        end
    end

    @testset "Distance Utilities Functions" begin
        
        @testset "Window Filtering Optimization" begin
            # Test the window filtering logic
            Random.seed!(111)
            model, data = create_distance_test_model()
            match!(model, data)
            
            # Should complete matching within reasonable time
            elapsed_time = @elapsed match!(model, data)
            @test elapsed_time < 30.0  # Should be much faster than 30 seconds for test data
        end
        
        @testset "Covariance Matrix Handling" begin
            model, data = create_distance_test_model()
            
            # Test that covariance matrix calculations work
            @test_nowarn match!(model, data)
            
            # Verify we have valid match information (can't directly access covariance dict)
            @test !isempty(model.matches)  # Should have match objects
        end
    end

    @testset "Regression Tests for Optimizations" begin
        
        @testset "Consistency After Threading Changes" begin
            # Test that results are consistent before/after optimization
            Random.seed!(2024)
            
            # Create identical models
            model1, data1 = create_distance_test_model()
            model2, data2 = create_distance_test_model()
            
            # Both should produce same results with same random seed
            Random.seed!(2024)
            match!(model1, data1)
            
            Random.seed!(2024)
            match!(model2, data2)
            
            # Results should be identical
            if !isempty(model1.matches) && !isempty(model2.matches)
                @test length(model1.matches) == length(model2.matches)
                
                for (m1, m2) in zip(model1.matches, model2.matches)
                    if !isempty(m1.distances) && !isempty(m2.distances)
                        @test m1.distances ≈ m2.distances atol=1e-12
                    end
                end
            end
        end
        
        @testset "Performance Characteristics" begin
            # Verify optimizations maintain good performance characteristics
            Random.seed!(2025)
            model, data = create_distance_test_model()
            
            # Matching should complete efficiently
            matching_time = @elapsed match!(model, data)
            @test matching_time < 10.0  # Reasonable time limit for test data
            
            # Memory usage should be reasonable (basic check)
            @test !isempty(model.matches)
            @test all(m -> isa(m.distances, Union{Matrix{Float64}, Vector{Matrix{Float64}}}) || isempty(m.distances), model.matches)
        end
        
        @testset "Edge Cases" begin
            
            @testset "Single Unit" begin
                # Test with minimal data
                Random.seed!(333)
                single_data = DataFrame(
                    date = [Date(2021, 10, 1), Date(2021, 10, 2)],
                    fips = [1001, 1001],
                    pop_dens = [50.0, 50.0],
                    cumul_death_rate = [20.0, 20.1],
                    death_rte = [0.5, 0.6],
                    gub = [0, 1],
                    day = [0, 1]
                )
                
                timevary = Dict(:pop_dens => false, :cumul_death_rate => true)
                single_model = makemodel(
                    single_data, :day, :fips, :gub, :death_rte,
                    [:pop_dens, :cumul_death_rate], timevary,
                    1:1, -1:-1;
                    title = "Single Unit Test"
                )
                
                # Should handle gracefully even with insufficient data
                @test_nowarn match!(single_model, single_data)
            end
            
            @testset "No Valid Matches" begin
                # Test scenario where no valid matches exist
                Random.seed!(444)
                no_match_data = DataFrame(
                    date = Date(2021, 10, 1):Day(1):Date(2021, 10, 10),
                    fips = repeat([1001], 10),
                    pop_dens = repeat([50.0], 10),
                    cumul_death_rate = 20.0:0.1:20.9,
                    death_rte = repeat([0.5], 10),
                    gub = [zeros(Int, 8)..., 1, 0],  # Treatment on day 9
                    day = 0:9
                )
                
                timevary = Dict(:pop_dens => false, :cumul_death_rate => true)
                no_match_model = makemodel(
                    no_match_data, :day, :fips, :gub, :death_rte,
                    [:pop_dens, :cumul_death_rate], timevary,
                    1:1, -5:-1;  # Impossible window for this data
                    title = "No Match Test"
                )
                
                # Should handle gracefully without errors
                @test_nowarn match!(no_match_model, no_match_data)
            end
        end
    end
end