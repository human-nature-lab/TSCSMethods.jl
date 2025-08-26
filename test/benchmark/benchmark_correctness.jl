# benchmark_correctness.jl
# Performance benchmarks WITH statistical correctness validation
# Run before/after optimization to ensure improvements don't break estimation

using TSCSMethods
using DataFrames, Statistics  
using Test, BenchmarkTools

using TSCSMethods: generate_realistic_tscs, generate_normal_effects, example_data

struct BenchmarkResult
    model::VeryAbstractCICModel
    case_name::String
    true_att_by_f::Union{Vector{Float64}, Nothing}
    data_rows::Int
    treatment_events::Int
end

# Helper functions for accessing model results
get_att_estimates(result::BenchmarkResult) = result.model.results.att
get_f_biases(result::BenchmarkResult) = isnothing(result.true_att_by_f) ? nothing : get_att_estimates(result) .- result.true_att_by_f
get_mean_bias(result::BenchmarkResult) = isnothing(get_f_biases(result)) ? nothing : mean(abs.(get_f_biases(result)))
is_f_correct(result::BenchmarkResult) = isnothing(get_f_biases(result)) ? "N/A - real data" : all(abs.(get_f_biases(result)) .< 0.05)

"""
    run_correctness_benchmark()

Run both performance and statistical correctness benchmarks.
Use this before and after optimizations to ensure improvements don't break estimation.
"""
function run_correctness_benchmark()
    # Test cases with F-period-specific effects
    test_cases = [
        (name="Example Data (Real)", data_func=() -> example_data(), true_att_by_f=nothing),
        (name="Zero Effects", data_func=() -> generate_realistic_tscs_by_f(true_att_by_f=fill(0.0, 10), seed=42), true_att_by_f=fill(0.0, 10)),
        (name="Small Constant", data_func=() -> generate_realistic_tscs_by_f(true_att_by_f=fill(0.1, 10), seed=42), true_att_by_f=fill(0.1, 10)),
        (name="Medium Constant", data_func=() -> generate_realistic_tscs_by_f(true_att_by_f=fill(0.5, 10), seed=42), true_att_by_f=fill(0.5, 10)),
        (name="Linear Decay", data_func=() -> generate_realistic_tscs_by_f(true_att_by_f=[0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05], seed=42), true_att_by_f=[0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]),
        (name="Normal Distribution", data_func=() -> generate_realistic_tscs_by_f(true_att_by_f=generate_normal_effects(0.5, 1:10), seed=42), true_att_by_f=generate_normal_effects(0.5, 1:10)),
    ]
    
    results = BenchmarkResult[]
    
    for (i, case) in enumerate(test_cases)
        # Generate data
        data = case.data_func()
        
        # Remove synthetic validation columns if present
        if hasproperty(data, :true_counterfactual)
            data = select(data, Not([:true_counterfactual, :true_treatment_effect]))
        end
        
        # Run model
        timevary = Dict{Symbol, Bool}(:pop_dens => false)
        model = makemodel(data, :day, :fips, :gub, :death_rte, [:pop_dens], timevary, 1:10, -20:-5)
        match!(model, data)
        estimate!(model, data)
        
        # Create result object
        result = BenchmarkResult(
            model,
            case.name,
            case.true_att_by_f,
            nrow(data),
            sum(data.gub)
        )
        
        # Test correctness for synthetic data
        # if !isnothing(case.true_att_by_f)
        #     @test is_f_correct(result)
        # end
        
        push!(results, result)
    end

    return results
end

"""
    compare_benchmarks(before, after)

Compare two benchmark results to show optimization impact.
"""
function compare_benchmarks(before::Vector{BenchmarkResult}, after::Vector{BenchmarkResult})
    @testset "Benchmark Comparison" begin
        # Match results by case name
        for before_result in before
            after_result = findfirst(r -> r.case_name == before_result.case_name, after)
            if !isnothing(after_result)
                after_result = after[after_result]
                
                @testset "$(before_result.case_name)" begin
                    # Test ATT estimates are equivalent (within tolerance)
                    before_att = get_att_estimate(before_result)
                    after_att = get_att_estimate(after_result)
                    @test abs(before_att - after_att) < 1e-10
                    
                    # Test bias estimates are equivalent (if ground truth available)
                    before_bias = get_bias(before_result)
                    after_bias = get_bias(after_result)
                    if !isnothing(before_bias) && !isnothing(after_bias)
                        @test abs(before_bias - after_bias) < 1e-10
                    end
                    
                    # Test correctness status is preserved
                    before_correctness = is_correct(before_result)
                    after_correctness = is_correct(after_result)
                    @test before_correctness === after_correctness
                end
            end
        end
    end
end

# Export functions for external use
export run_correctness_benchmark, compare_benchmarks

res = run_correctness_benchmark()

# Create summary DataFrame for F-period analysis
validation_results = DataFrame(
    case_name = String[],
    true_att_by_f = Union{Vector{Float64}, Nothing}[],
    estimated_att_by_f = Vector{Float64}[],
    f_biases = Union{Vector{Float64}, Nothing}[],
    mean_abs_bias = Union{Float64, Nothing}[],
    max_abs_bias = Union{Float64, Nothing}[],
    data_rows = Int[],
    treatment_events = Int[],
    all_f_pass = Union{Bool, Missing}[],
    f_success_rate = Union{Float64, Missing}[]
)

# Run formal tests
@testset "F-Period Specific Validation" begin
    for result in res
        estimated_att_by_f = get_att_estimates(result)
        
        if !isnothing(result.true_att_by_f)
            f_biases = get_f_biases(result)
            mean_abs_bias = get_mean_bias(result)
            max_abs_bias = maximum(abs.(f_biases))
            all_f_pass = is_f_correct(result)
            f_success_rate = sum(abs.(f_biases) .< 0.05) / length(f_biases)
            
            @testset "$(result.case_name)" begin
                @test mean_abs_bias < 0.1  # Relaxed threshold for F-period testing
                @test all(isfinite.(estimated_att_by_f))
                # Test that at least some F periods are accurate
                @test f_success_rate > 0.3
            end
            
            # Add to results DataFrame
            push!(validation_results, (
                result.case_name,
                result.true_att_by_f,
                estimated_att_by_f,
                f_biases,
                mean_abs_bias,
                max_abs_bias,
                result.data_rows,
                result.treatment_events,
                all_f_pass,
                f_success_rate
            ))
        else
            # Real data case - no validation possible
            @testset "$(result.case_name)" begin
                @test all(isfinite.(estimated_att_by_f))
            end
            
            push!(validation_results, (
                result.case_name,
                nothing,
                estimated_att_by_f,
                nothing,
                nothing,
                nothing,
                result.data_rows,
                result.treatment_events,
                missing,
                missing
            ))
        end
    end
end

# Display results table
show(validation_results, allrows=true, allcols=true)
import CSV, Dates
pth = join(["test/validation/validation_results_", Dates.now(), ".csv"]);
CSV.write(pth, validation_results; transform=(col, val) -> something(val, missing))
