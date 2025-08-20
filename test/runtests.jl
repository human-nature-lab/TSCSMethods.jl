using TSCSMethods
using Test
using DataFrames
using Dates
using Accessors

@testset "TSCSMethods.jl" begin
    include("test_workflow.jl")  # Main workflow tests first
    include("test_simple.jl")
    include("test_construction.jl")
    include("test_matching.jl")
    include("test_balancing.jl")
    include("test_estimation.jl")
    include("test_utilities.jl")
end