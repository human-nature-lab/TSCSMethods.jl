# benchmark_correctness.jl
# Performance benchmarks WITH statistical correctness validation
# Run before/after optimization to ensure improvements don't break estimation

using TSCSMethods
using DataFrames, Statistics  
using BenchmarkTools
include("data_generation.jl")

"""
    run_correctness_benchmark()

Run both performance and statistical correctness benchmarks.
Use this before and after optimizations to ensure improvements don't break estimation.
"""
function run_correctness_benchmark()
    println("🚀 PERFORMANCE + CORRECTNESS BENCHMARK")
    println("=" ^ 60)
    
    # Test cases
    test_cases = [
        (name="Example Data (Real)", data_func=() -> example_data(), true_att=nothing),
        (name="Realistic Synthetic (Zero)", data_func=() -> generate_realistic_tscs(true_att=0.0, seed=42), true_att=0.0),
        (name="Realistic Synthetic (Negative)", data_func=() -> generate_realistic_tscs(true_att=-0.02, seed=123), true_att=-0.02),
    ]
    
    results = Dict{String, Any}()
    
    for (i, case) in enumerate(test_cases)
        println("\n📊 Test Case $i: $(case.name)")
        println("-" ^ 40)
        
        # Generate data
        data = case.data_func()
        
        # Remove synthetic validation columns if present
        if hasproperty(data, :true_counterfactual)
            data = select(data, Not([:true_counterfactual, :true_treatment_effect]))
        end
        
        println("   Data: $(nrow(data)) rows, $(sum(data.gub)) treatment events")
        
        # Performance benchmark
        timevary = Dict{Symbol, Bool}(:pop_dens => false)
        
        benchmark_result = @benchmark begin
            model = makemodel($data, :day, :fips, :gub, :death_rte, [:pop_dens], $timevary, 1:10, -20:-5)
            match!(model, $data)
            estimate!(model, $data)
            model.results.att[1]  # Return first ATT estimate
        end seconds=30 samples=10
        
        # Extract performance metrics
        median_time = median(benchmark_result.times) / 1e9  # Convert to seconds
        memory_mb = median(benchmark_result.memory) / 1024^2  # Convert to MB
        
        println("   ⏱️  Time: $(round(median_time, digits=3))s")
        println("   💾 Memory: $(round(memory_mb, digits=1)) MB")
        
        # Statistical correctness check
        model = makemodel(data, :day, :fips, :gub, :death_rte, [:pop_dens], timevary, 1:10, -20:-5)
        match!(model, data)
        estimate!(model, data)
        
        att_estimate = model.results.att[1]
        
        if case.true_att !== nothing
            bias = att_estimate - case.true_att
            println("   🎯 ATT: $(round(att_estimate, digits=4)) (true: $(case.true_att), bias: $(round(bias, digits=4)))")
            
            # Correctness validation
            if abs(bias) < 0.05
                println("   ✅ Statistical correctness: PASS")
                correctness_status = "PASS"
            else
                println("   ❌ Statistical correctness: FAIL (bias too large)")
                correctness_status = "FAIL"
            end
        else
            println("   📈 ATT: $(round(att_estimate, digits=4)) (real data - no ground truth)")
            correctness_status = "N/A"
        end
        
        # Store results
        results[case.name] = Dict(
            "median_time_s" => median_time,
            "memory_mb" => memory_mb,
            "att_estimate" => att_estimate,
            "true_att" => case.true_att,
            "bias" => case.true_att !== nothing ? att_estimate - case.true_att : nothing,
            "correctness_status" => correctness_status,
            "data_rows" => nrow(data),
            "treatment_events" => sum(data.gub)
        )
    end
    
    # Summary
    println("\n📋 BENCHMARK SUMMARY")
    println("=" ^ 60)
    
    total_time = sum(r["median_time_s"] for r in values(results))
    max_memory = maximum(r["memory_mb"] for r in values(results))
    correctness_passes = sum(r["correctness_status"] == "PASS" for r in values(results))
    total_correctness_tests = sum(r["correctness_status"] != "N/A" for r in values(results))
    
    println("⏱️  Total time: $(round(total_time, digits=3))s")
    println("💾 Peak memory: $(round(max_memory, digits=1)) MB") 
    println("✅ Correctness: $correctness_passes/$total_correctness_tests tests passed")
    
    if correctness_passes == total_correctness_tests
        println("🎉 All correctness tests passed - safe to optimize!")
    else
        println("⚠️  Some correctness tests failed - investigate before optimizing")
    end
    
    return results
end

"""
    compare_benchmarks(before, after)

Compare two benchmark results to show optimization impact.
"""
function compare_benchmarks(before::Dict, after::Dict)
    println("\n📊 OPTIMIZATION IMPACT ANALYSIS")
    println("=" ^ 60)
    
    for test_name in keys(before)
        if test_name in keys(after)
            println("\n🧪 $test_name:")
            
            before_time = before[test_name]["median_time_s"]
            after_time = after[test_name]["median_time_s"]
            time_improvement = (before_time - after_time) / before_time * 100
            
            before_memory = before[test_name]["memory_mb"] 
            after_memory = after[test_name]["memory_mb"]
            memory_improvement = (before_memory - after_memory) / before_memory * 100
            
            println("   ⏱️  Time: $(round(before_time, digits=3))s → $(round(after_time, digits=3))s ($(round(time_improvement, digits=1))% improvement)")
            println("   💾 Memory: $(round(before_memory, digits=1))MB → $(round(after_memory, digits=1))MB ($(round(memory_improvement, digits=1))% improvement)")
            
            # Check correctness preservation
            before_correctness = before[test_name]["correctness_status"]
            after_correctness = after[test_name]["correctness_status"] 
            
            if before_correctness == after_correctness == "PASS"
                println("   ✅ Correctness: Maintained")
            elseif before_correctness == "PASS" && after_correctness == "FAIL"
                println("   ❌ Correctness: BROKEN by optimization!")
            elseif before_correctness == "FAIL" && after_correctness == "PASS"
                println("   ✅ Correctness: Fixed by optimization")
            else
                println("   📝 Correctness: $(before_correctness) → $(after_correctness)")
            end
        end
    end
end

# Export functions for external use
export run_correctness_benchmark, compare_benchmarks