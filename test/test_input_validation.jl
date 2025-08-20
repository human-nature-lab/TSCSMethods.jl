"""
Test suite for input validation and error handling in distance calculation functions.

This test suite ensures that all distance calculation functions properly validate their inputs
and provide meaningful error messages for invalid data, while maintaining correctness
for valid inputs.
"""

using Test
using TSCSMethods
using Statistics
using LinearAlgebra

@testset "Input Validation & Error Handling Tests" begin

    @testset "get_thread_storage validation" begin
        @testset "Invalid dimensions" begin
            @test_throws ArgumentError TSCSMethods.get_thread_storage(0, 5)
            @test_throws ArgumentError TSCSMethods.get_thread_storage(-1, 5)
            @test_throws ArgumentError TSCSMethods.get_thread_storage(10, 0)
            @test_throws ArgumentError TSCSMethods.get_thread_storage(10, -1)
        end
        
        @testset "Valid inputs" begin
            storage = TSCSMethods.get_thread_storage(10, 3)
            @test storage isa TSCSMethods.ThreadLocalDistanceStorage
            @test length(storage.dtots_float) == 4  # 3 covariates + 1 for mahalanobis
            @test length(storage.dtots_missing) == 4
            @test storage.max_times >= 10
        end
        
        @testset "Large dimension warning" begin
            # Capture warnings
            @test_logs (:warn, r"Large n_times") TSCSMethods.get_thread_storage(15000, 3)
        end
    end

    @testset "samplecovar validation" begin
        @testset "Empty inputs" begin
            @test_throws ArgumentError TSCSMethods.samplecovar(Int[], Matrix{Float64}(undef, 0, 2))
            @test_throws ArgumentError TSCSMethods.samplecovar([1, 2], Matrix{Float64}(undef, 2, 0))
        end
        
        @testset "Dimension mismatch" begin
            @test_throws DimensionMismatch TSCSMethods.samplecovar([1, 2, 3], randn(2, 3))
            @test_throws DimensionMismatch TSCSMethods.samplecovar([1, 2], randn(3, 2))
        end
        
        @testset "Valid inputs" begin
            dat_t = [1, 1, 2, 2, 3, 3]
            cdat = randn(6, 3)
            result = TSCSMethods.samplecovar(dat_t, cdat)
            @test result isa Dict{Int64, Matrix{Float64}}
            @test length(result) == 3  # Three time periods
        end
    end

    @testset "alldistances! validation" begin
        @testset "Empty inputs" begin
            @test_throws ArgumentError TSCSMethods.alldistances!(
                Vector{Vector{Float64}}(), Dict{Int, Matrix{Float64}}(), 
                Vector{Float64}[], Vector{Float64}[], Int[]
            )
            
            dtotals = [Vector{Float64}(undef, 0)]
            @test_throws ArgumentError TSCSMethods.alldistances!(
                dtotals, Dict{Int, Matrix{Float64}}(), 
                Vector{Float64}[], Vector{Float64}[], Int[]
            )
        end
        
        @testset "Dimension mismatches" begin
            dtotals = [Vector{Float64}(undef, 3), Vector{Float64}(undef, 3)]
            Σinvdict = Dict(1 => Matrix{Float64}(I, 2, 2))
            
            # Length mismatches between xrows, yrows, γtimes
            @test_throws DimensionMismatch TSCSMethods.alldistances!(
                dtotals, Σinvdict,
                [randn(2), randn(2)],  # 2 time points
                [randn(2)],           # 1 time point - mismatch!
                [1, 2]                # 2 time points
            )
            
            # dtotals length mismatch with γtimes
            @test_throws DimensionMismatch TSCSMethods.alldistances!(
                [Vector{Float64}(undef, 2)],  # Wrong length
                Σinvdict,
                [randn(2), randn(2), randn(2)],  # 3 time points
                [randn(2), randn(2), randn(2)],
                [1, 2, 3]
            )
        end
        
        @testset "Valid computation" begin
            # Setup valid test data
            n_times = 5
            n_covs = 2
            dtotals = [Vector{Float64}(undef, n_times) for _ in 1:(n_covs+1)]
            
            Σinvdict = Dict(
                i => Matrix{Float64}(I, n_covs, n_covs) for i in 1:n_times
            )
            
            xrows = [randn(n_covs) for _ in 1:n_times]
            yrows = [randn(n_covs) for _ in 1:n_times]
            γtimes = 1:n_times
            
            # Should not throw
            @test_nowarn TSCSMethods.alldistances!(dtotals, Σinvdict, xrows, yrows, γtimes)
            
            # Check that distances were computed
            @test all(isfinite, dtotals[1])  # Mahalanobis distances
            @test all(isfinite, dtotals[2])  # Covariate 1 distances
            @test all(isfinite, dtotals[3])  # Covariate 2 distances
        end
    end

    @testset "distaveraging! validation (sliding window)" begin
        @testset "Empty inputs" begin
            distances = Matrix{Float64}[]
            dtots = Vector{Vector{Float64}}()
            accums = Int[]
            
            @test_throws ArgumentError TSCSMethods.distaveraging!(
                distances, dtots, accums, Int[], 1:5, 1, 1
            )
        end
        
        @testset "Dimension mismatches" begin
            n_times = 5
            distances = [Matrix{Float64}(undef, 2, 2), Matrix{Float64}(undef, 2, 2)]
            dtots = [Vector{Float64}(undef, n_times), Vector{Float64}(undef, 3)]  # Mismatch!
            accums = [0, 0]
            γtimes = 1:n_times
            
            @test_throws DimensionMismatch TSCSMethods.distaveraging!(
                distances, dtots, accums, γtimes, 1:3, 1, 1
            )
        end
        
        @testset "Invalid indices" begin
            distances = [Matrix{Float64}(undef, 2, 2)]
            dtots = [Vector{Float64}(undef, 5)]
            accums = [0]
            
            @test_throws BoundsError TSCSMethods.distaveraging!(
                distances, dtots, accums, 1:5, 1:3, 0, 1  # φ = 0 invalid
            )
            
            @test_throws BoundsError TSCSMethods.distaveraging!(
                distances, dtots, accums, 1:5, 1:3, 1, 0  # m = 0 invalid
            )
            
            @test_throws BoundsError TSCSMethods.distaveraging!(
                distances, dtots, accums, 1:5, 1:3, 3, 1  # φ = 3 > matrix size
            )
        end
        
        @testset "Valid computation" begin
            n_times = 10
            distances = [Matrix{Float64}(undef, 3, 3) for _ in 1:2]
            dtots = [randn(n_times), randn(n_times)]
            accums = [0, 0]
            γtimes = 1:n_times
            fw = 3:7  # Window from time 3 to 7
            
            @test_nowarn TSCSMethods.distaveraging!(distances, dtots, accums, γtimes, fw, 1, 1)
            
            # Check results are finite and reasonable
            @test isfinite(distances[1][1, 1])
            @test isfinite(distances[2][1, 1])
        end
    end

    @testset "distaveraging! validation (fixed window)" begin
        @testset "Empty window" begin
            drow = [0.0, 0.0]
            dtots = [Vector{Float64}(undef, 5), Vector{Float64}(undef, 5)]
            accums = [0, 0]
            
            @test_throws ArgumentError TSCSMethods.distaveraging!(
                drow, dtots, accums, 1:5, Int[]  # Empty fw
            )
        end
        
        @testset "Valid computation with missing data" begin
            n_times = 8
            drow = [0.0, 0.0]
            dtots = [
                Vector{Union{Float64, Missing}}([1.0, 2.0, missing, 4.0, 5.0, 6.0, 7.0, 8.0]),
                Vector{Union{Float64, Missing}}([2.0, missing, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0])
            ]
            accums = [0, 0]
            γtimes = 1:n_times
            fw = 2:6  # Window from time 2 to 6
            
            @test_nowarn TSCSMethods.distaveraging!(drow, dtots, accums, γtimes, fw)
            
            # Should compute averages excluding missing values
            @test isfinite(drow[1])
            @test isfinite(drow[2]) 
            @test drow[1] ≈ (2.0 + 4.0 + 5.0 + 6.0) / 4  # Excluding missing at index 3
            @test drow[2] ≈ (6.0 + 8.0 + 10.0 + 12.0) / 4  # Excluding missing at index 2
        end
    end

    @testset "mahaveraging validation" begin
        @testset "Dimension mismatches" begin
            @test_throws DimensionMismatch TSCSMethods.mahaveraging(
                [1.0, 2.0], [1, 2, 3], 1:2  # Length mismatch
            )
        end
        
        @testset "Empty window" begin
            @test_throws ArgumentError TSCSMethods.mahaveraging(
                [1.0, 2.0], [1, 2], Int[]
            )
        end
        
        @testset "No valid observations warning" begin
            # All times outside window
            @test_logs (:warn, r"No valid observations") TSCSMethods.mahaveraging(
                [1.0, 2.0], [10, 11], 1:2
            )
        end
        
        @testset "Valid computation" begin
            mahas = [1.0, 2.0, 3.0, 4.0]
            T = [1, 2, 3, 4]
            fw = 2:3
            
            result = TSCSMethods.mahaveraging(mahas, T, fw)
            @test result ≈ (2.0 + 3.0) / 2
        end
        
        @testset "Missing data handling" begin
            mahas = Union{Float64, Missing}[1.0, missing, 3.0, 4.0]
            T = [1, 2, 3, 4]
            fw = 1:4
            
            result = TSCSMethods.mahaveraging(mahas, T, fw)
            @test result ≈ (1.0 + 3.0 + 4.0) / 3  # Excluding missing
        end
    end


    @testset "distances_calculate! validation" begin
        @testset "Empty observations" begin
            matches = []
            observations = []
            
            @test_logs (:warn, "No observations to process") TSCSMethods.distances_calculate!(
                matches, observations, Int[], String[], nothing, nothing, 
                1, 1, 2, Dict{Int, Matrix{Float64}}()
            )
        end
        
        @testset "Length mismatches" begin
            # Create minimal valid structures for testing
            matches = [TSCSMethods.Tob(mus=Matrix{Bool}(undef, 0, 0), distances=Matrix{Float64}(undef, 0, 0), ranks=Dict{Int, Vector{Int}}())]
            observations = [1, 2]  # Length mismatch
            
            @test_throws ArgumentError TSCSMethods.distances_calculate!(
                matches, observations, [1], ["var1"], nothing, nothing,
                1, 1, 2, Dict(1 => Matrix{Float64}(I, 1, 1))
            )
        end
        
        @testset "Invalid time bounds" begin
            matches = [TSCSMethods.Tob(mus=Matrix{Bool}(undef, 0, 0), distances=Matrix{Float64}(undef, 0, 0), ranks=Dict{Int, Vector{Int}}())]
            observations = [1]
            
            @test_throws ArgumentError TSCSMethods.distances_calculate!(
                matches, observations, [1], ["var1"], nothing, nothing,
                1, 3, 2, Dict(1 => Matrix{Float64}(I, 1, 1))  # Lmin=3 > Lmax=2
            )
        end
        
        @testset "Sliding window not implemented" begin
            matches = [TSCSMethods.Tob(mus=Matrix{Bool}(undef, 0, 0), distances=Matrix{Float64}(undef, 0, 0), ranks=Dict{Int, Vector{Int}}())]
            observations = [1]
            
            @test_throws ErrorException TSCSMethods.distances_calculate!(
                matches, observations, [1], ["var1"], nothing, nothing,
                1, 1, 2, Dict(1 => Matrix{Float64}(I, 1, 1)); sliding = true
            )
        end
    end

    @testset "Edge Cases and Robustness" begin
        @testset "Single observation" begin
            # Test that functions handle single-element arrays correctly
            @test_nowarn TSCSMethods.mahaveraging([1.0], [1], 1:1)
        end
        
        @testset "Large inputs" begin
            # Test with reasonably large inputs to ensure no overflow
            large_n = 1000
            mahas = randn(large_n)
            T = 1:large_n
            fw = 100:900
            
            result = TSCSMethods.mahaveraging(mahas, T, fw)
            @test isfinite(result)
        end
        
        @testset "Extreme values" begin
            # Test with very large/small finite values
            mahas = [1e-10, 1e10, -1e10, 1e-10]
            T = [1, 2, 3, 4]
            fw = 1:4
            
            result = TSCSMethods.mahaveraging(mahas, T, fw)
            @test isfinite(result)
        end
    end
end