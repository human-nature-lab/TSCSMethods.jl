using TSCSMethods
using Test
using DataFrames
include(joinpath(@__DIR__, "..", "support", "simulate_tscs.jl"))

"""
Compute naive event-time DiD ATT under randomized assumptions (will be biased under confounding).
"""
function naive_event_time_att(df::DataFrame, treated_ids, t0::Dict{Int,Int}, F)
    never = setdiff(unique(df.id), treated_ids)
    att = zeros(Float64, length(F))
    for (k, f) in enumerate(F)
        s = 0.0; c = 0
        for u in treated_ids
            tt = t0[u]
            y_t = df[(df.id .== u) .& (df.t .== tt + f), :y][1]
            y_r = df[(df.id .== u) .& (df.t .== tt - 1), :y][1]
            yc_t = mean(df[(in.(df.id, Ref(never))) .& (df.t .== tt + f), :y])
            yc_r = mean(df[(in.(df.id, Ref(never))) .& (df.t .== tt - 1), :y])
            s += (y_t - y_r) - (yc_t - yc_r)
            c += 1
        end
        att[k] = s / max(c, 1)
    end
    return att
end

@testset "Confounded DGP bias reduction" begin
    F = 1:5    # Post-treatment periods (positive)
    L = -5:-1  # Pre-treatment periods (negative)  
    delta = [0.01, 0.02, 0.05, 0.03, 0.01]  # Small true effects for bias detection

    df, treated_ids, t0 = simulate_confounded(
        N=600, T=200, treated_share=0.10, F=F, L=L, delta=delta,
        phi=0.0, rho=0.5, sigma_y=0.2, beta1=0.3, beta2=0.5, seed=321,
        conf_strength=0.7, pretrend_strength=0.02
    )

    naive = naive_event_time_att(df, treated_ids, t0, F)
    naive_bias = abs.(naive .- delta)
    naive_mae = mean(naive_bias)

    timevary = Dict(:x1 => false, :x2 => true)
    model = makemodel(df, :t, :id, :gub, :y, [:x1, :x2], timevary, F, L)
    model.iterations = 100
    match!(model, df)
    # attempt caliper/refinement based balancing to reduce selection bias
    try
        autobalance(model, df; threshold=0.1, refinementnum=3, calmin=0.1, step=0.05, doestimate=false, verbose=false)
    catch
        # proceed if autobalance not applicable
    end
    estimate!(model, df; dobayesfactor=false)
    est = collect(model.results.att)
    match_bias = abs.(est .- delta)
    match_mae = mean(match_bias)

    @info "Bias reduction (naive vs matched)" naive_mae match_mae naive est delta

    # PR-safe thresholds per specification
    bias_reduction = (naive_mae - match_mae) / naive_mae
    @test bias_reduction >= 0.30 "Mean absolute bias reduced by ≥30% vs naive"
    @test match_mae <= 0.05 "Post-matching mean absolute bias ≤ 0.05"
end

