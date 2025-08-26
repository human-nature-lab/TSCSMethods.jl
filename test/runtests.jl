using TSCSMethods
using Test
using DataFrames
using Dates
using Accessors
using Random
using StatsBase

@testset "Backend tests" begin
    include("test_simple.jl")
    include("test_construction.jl")
    include("test_matching.jl")
    include("test_distance.jl")  # Distance calculation optimization tests
    include("test_balancing.jl")
    include("test_estimation.jl")
    include("test_utilities.jl")
end

@testset "Correctness tests" begin
    include("test_noiseless_exact.jl")
    include("test_synthetic_known_effects.jl")
    include("test_confounded_dgp.jl")
    include("test_coverage_type1.jl")
end
