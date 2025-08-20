#!/usr/bin/env julia

# Simple compatibility test for TSCSMethods.jl
println("Testing TSCSMethods.jl compatibility...")
println("Julia version: ", VERSION)

using Pkg
Pkg.activate(".")

try
    println("Loading TSCSMethods...")
    using TSCSMethods, DataFrames
    
    println("Testing core workflow...")
    # Quick test of core functionality
    dat = example_data(n_units=5, n_days=30, seed=42)
    model = makemodel(dat, :day, :fips, :gub, :death_rte, [:pop_dens], Dict(:pop_dens => false), 5:8, -10:-5)
    match!(model, dat)
    balance!(model, dat)
    estimate!(model, dat; dobayesfactor=false)
    
    println("✅ Core workflow successful on Julia $(VERSION)")
    println("✅ Generated $(nrow(model.results)) ATT estimates")
    
catch e
    println("❌ Compatibility test failed: ", e)
    exit(1)
end

println("🎉 TSCSMethods.jl compatible with Julia $(VERSION)")