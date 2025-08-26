using TSCSMethods
using Test
using DataFrames
using Dates
using Accessors
using Random
using StatsBase

@testset "Backend tests" begin
    include("integration/test_simple.jl")
    include("unit/test_construction.jl")
    include("unit/test_matching.jl")
    include("unit/test_distance.jl")  # Distance calculation optimization tests
    include("unit/test_balancing.jl")
    include("unit/test_estimation.jl")
    include("unit/test_utilities.jl")
end

@testset "Correctness tests" begin
    include("correctness/test_noiseless_exact.jl")
    include("correctness/test_synthetic_known_effects.jl")
    include("correctness/test_confounded_dgp.jl")
    include("correctness/test_coverage_type1.jl")
end
