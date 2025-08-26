using TSCSMethods
using Test
using DataFrames
using Random
using Statistics
using StatsBase: sample

# Simple TSCS simulator with known per-f effects
function simulate_simple(; N::Int=600, T::Int=200, treated_share::Float64=0.05,
    F::UnitRange{Int}=1:8, L::UnitRange{Int}=-30:-1, delta::Vector{Float64}=[0.00, 0.02, 0.04, 0.05, 0.05, 0.03, 0.02, 0.00],
    φ::Float64=0.0, ρ::Float64=0.5, σy::Float64=0.2, β1::Float64=0.3, β2::Float64=0.5, seed::Int=42)

    @assert length(delta) == length(F) "delta length must match |F|"
    Random.seed!(seed)

    ids = collect(1:N)
    tgrid = collect(1:T)
    n_treated = max(1, round(Int, treated_share*N))
    treated_ids = sort(sample(ids, n_treated; replace=false))

    # Ensure event day allows full [t0+Lmin, t0+fmax) window within [1, T]
    Lmin = minimum(L); fmax = maximum(F)
    # Ensure: t0 + Lmin >= 1  and  t0 + fmax <= T
    t0_min = 1 - Lmin         # lower bound for valid pre-periods
    t0_max = T - fmax         # ensure t0 + fmax <= T
    @assert t0_min < t0_max "Time grid too short for chosen L and F"
    t0_map = Dict(u => rand(t0_min:t0_max) for u in treated_ids)

    # Covariates
    x1 = Dict(u => randn() for u in ids) # time-invariant
    x2 = Dict{Tuple{Int,Int},Float64}()
    for u in ids
        prev = randn()
        for t in tgrid
            prev = ρ*prev + randn()*0.5
            x2[u,t] = prev
        end
    end

    # Unit/time effects
    α = Dict(u => randn()*0.5 for u in ids)
    τ = [randn()*0.2 for _ in tgrid]

    # Outcome baseline AR(1)
    y = Dict{Tuple{Int,Int},Float64}()
    for u in ids
        y_prev = randn()*σy
        for (ti, t) in enumerate(tgrid)
            μ = α[u] + τ[ti] + β1*x1[u] + β2*x2[u,t]
            y_curr = μ + φ*y_prev + randn()*σy
            y[u,t] = y_curr
            y_prev = y_curr
        end
    end

    # Apply known effects at event-time offsets
    for u in treated_ids
        t0 = t0_map[u]
        for (j, f) in enumerate(F)
            tt = t0 + f
            y[u, tt] = y[u, tt] + delta[j]
        end
    end

    # Build DataFrame
    rows = Vector{NamedTuple}(undef, N*T)
    k = 0
    for u in ids
        t0 = get(t0_map, u, -10^9) # non-treated units have no event day
        for t in tgrid
            k += 1
            rows[k] = (id=u, t=t, gub=Int(t == t0), y=y[u,t], x1=x1[u], x2=x2[u,t])
        end
    end
    return DataFrame(rows), treated_ids, t0_map
end

@testset "Synthetic data recovers known per-f ATTs" begin
    # Ground-truth effects for F
    F = 1:8
    L = -30:-1
    delta = [0.00, 0.02, 0.04, 0.05, 0.05, 0.03, 0.02, 0.00]

    df, treated_ids, t0 = simulate_simple(; N=1000, T=220, treated_share=0.10, σy=0.15, F=F, L=L, delta=delta, seed=123)

    # Model setup
    timevary = Dict(:x1 => false, :x2 => true)
    model = makemodel(df, :t, :id, :gub, :y, [:x1, :x2], timevary, F, L)

    # Matching and estimation (no balancing needed; randomized DGP)
    match!(model, df)
    model.iterations = 100  # Runtime control as per spec
    estimate!(model, df; dobayesfactor=false)

    est = copy(model.results.att)
    # Inspect support (treated counts and matches) for context on stability
    supports = (treated = copy(model.results.treated), matches = copy(model.results.matches))

    # Naive oracle: event-time DiD under randomization to validate DGP
    naive = fill(0.0, length(F))
    never = setdiff(unique(df.id), treated_ids)
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
        naive[k] = s / max(c, 1)
    end

    # Compare to ground truth with reasonable tolerances
    mae = mean(abs.(est .- delta))
    mxe = maximum(abs.(est .- delta))

    @info "Known delta vs estimated ATT" delta est mae mxe supports naive

    @test mae ≤ 0.025
    @test mxe ≤ 0.06
end
