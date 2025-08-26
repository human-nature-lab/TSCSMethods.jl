using TSCSMethods
using Test
using DataFrames

include(joinpath(@__DIR__, "..", "support", "simulate_tscs.jl"))

@testset "Noiseless exact recovery" begin
    # Choose modest windows with valid support
    F = 1:6
    L = -15:-1
    delta = [0.00, 0.02, 0.04, 0.05, 0.03, 0.00]

    # Generate noiseless data: sigma_y=0, phi=0; set beta2=0 to remove time-varying covariate effect
    df, treated_ids, t0 = simulate_noiseless(
        N=400, T=160, treated_share=0.10,
        F=F, L=L, delta=delta, rho=0.5, beta1=0.4, beta2=0.0, seed=777
    )

    # Oracle ATT equals injected deltas exactly in noiseless mode
    true_att = oracle_att(df, treated_ids, t0, F)
    @test maximum(abs.(true_att .- delta)) ≤ 1e-12

    # Estimation pipeline should recover ATT exactly (up to machine precision)
    timevary = Dict(:x1 => false, :x2 => true)
    model = makemodel(df, :t, :id, :gub, :y, [:x1, :x2], timevary, F, L)
    match!(model, df)
    estimate!(model, df; dobayesfactor=false)

    est = collect(model.results.att)
    @test maximum(abs.(est .- delta)) ≤ 1e-3
end

