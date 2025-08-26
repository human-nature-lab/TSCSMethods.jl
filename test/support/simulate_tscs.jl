using DataFrames
using Random
using Statistics
using StatsBase: sample, Weights

"""
    simulate_tscs(; N, T, treated_share, F, L, delta, phi, rho, sigma_y, beta1, beta2, seed,
                    confounded::Bool=false, conf_strength::Float64=0.7,
                    pretrend::Bool=false, pretrend_strength::Float64=0.02)

Generate TSCS data with single-day event-time treatment and known per-f effects.

Returns: (df::DataFrame, treated_ids::Vector{Int}, t0_map::Dict{Int,Int}) with columns
`id, t, gub, y, x1, x2`.
"""
function simulate_tscs(; N::Int, T::Int, treated_share::Float64, F::UnitRange{Int}, L::UnitRange{Int},
    delta::Vector{Float64}, phi::Float64, rho::Float64, sigma_y::Float64, beta1::Float64, beta2::Float64,
    seed::Int, confounded::Bool=false, conf_strength::Float64=0.7,
    pretrend::Bool=false, pretrend_strength::Float64=0.02,
    noiseless::Bool=false)

    @assert length(delta) == length(F) "delta length must match |F|"
    Random.seed!(seed)

    ids = collect(1:N)
    tgrid = collect(1:T)

    # Covariates
    x1 = Dict(u => randn() for u in ids) # time-invariant
    x2 = Dict{Tuple{Int,Int},Float64}()
    for u in ids
        prev = randn()
        for t in tgrid
            prev = rho*prev + randn()*0.5
            x2[u,t] = prev
        end
    end

    # Unit/time effects
    α = Dict(u => randn()*0.5 for u in ids)
    τ = [randn()*0.2 for _ in tgrid]

    # Choose treated units
    n_treated = max(1, round(Int, treated_share*N))
    treated_ids = if confounded
        # Probability weight increases with x1 mean and x2 level at middle time
        w = [exp(conf_strength * (0.7*x1[u])) for u in ids]
        sample(ids, Weights(w), n_treated; replace=false)
    else
        sample(ids, n_treated; replace=false)
    end |> sort

    # Event times with window safety
    Lmin = minimum(L); fmax = maximum(F)
    t0_min = 1 - Lmin
    t0_max = T - fmax
    @assert t0_min < t0_max "Time grid too short for chosen L and F"

    if confounded
        # weak dependence of timing on x1
        probs = [exp(0.1*x1[u]) for u in treated_ids]
        t0_map = Dict{Int,Int}()
        for (u, p) in zip(treated_ids, probs)
            # sample around the center with slight shift
            base = rand(t0_min:t0_max)
            t0_map[u] = clamp(base, t0_min, t0_max)
        end
    else
        t0_map = Dict(u => rand(t0_min:t0_max) for u in treated_ids)
    end

    # Pretrend slopes (only affect t < t0)
    slope = Dict(u => (pretrend ? randn()*pretrend_strength : 0.0) for u in ids)

    # Outcome baseline AR(1)
    y = Dict{Tuple{Int,Int},Float64}()
    y0 = Dict{Tuple{Int,Int},Float64}()  # noise-inclusive counterfactual (no treatment effect)
    # Apply noiseless overrides
    σy = noiseless ? 0.0 : sigma_y
    ϕ = noiseless ? 0.0 : phi
    for u in ids
        y_prev = randn()*σy
        for (ti, t) in enumerate(tgrid)
            μ = α[u] + τ[ti] + beta1*x1[u] + beta2*x2[u,t]
            y_curr = μ + ϕ*y_prev + randn()*σy
            # add pretrend before event time if defined
            if pretrend
                t0u = get(t0_map, u, T+10)
                if t < t0u
                    y_curr += slope[u] * (t - t0u)
                end
            end
            y[u,t] = y_curr
            y0[u,t] = y_curr
            y_prev = y_curr
        end
    end

    # Apply known effects at event-time offsets
    for u in treated_ids
        t0 = t0_map[u]
        for (j, f) in enumerate(F)
            tt = t0 + f
            if 1 <= tt <= T  # Ensure tt is within valid time range
                y[u, tt] = y[u, tt] + delta[j]
            end
        end
    end

    # Build DataFrame
    rows = Vector{NamedTuple}(undef, N*T)
    k = 0
    for u in ids
        t0 = get(t0_map, u, -10^9)
        for t in tgrid
            k += 1
            rows[k] = (id=u, t=t, gub=Int(t == t0), y=y[u,t], y0=y0[u,t], x1=x1[u], x2=x2[u,t])
        end
    end
    return DataFrame(rows), treated_ids, t0_map
end

simulate_randomized(; kwargs...) = simulate_tscs(; confounded=false, pretrend=false, kwargs...)
simulate_confounded(; kwargs...) = simulate_tscs(; confounded=true, pretrend=true, kwargs...)
simulate_null(; kwargs...) = simulate_tscs(; kwargs..., delta=zeros(length(kwargs[:F])))

"""
    simulate_noiseless(; kwargs...)

Convenience wrapper for exact-recovery sanity checks. Sets `sigma_y=0` and `phi=0` to remove stochastic noise.
"""
simulate_noiseless(; kwargs...) = simulate_tscs(; kwargs..., sigma_y=0.0, phi=0.0, noiseless=true)

"""
    oracle_att(df, treated_ids, t0_map, F) -> Vector{Float64}

Compute the oracle ATT using the simulator's noise-inclusive counterfactual `y0`.
Returns mean over treated events of (y − y0) at each event-time offset f ∈ F.
Requires that `df` contain columns `:y` and `:y0` (provided by simulate_tscs).
"""
function oracle_att(df::DataFrame, treated_ids, t0_map::Dict{Int,Int}, F::UnitRange{Int})
    att = zeros(Float64, length(F))
    counts = zeros(Int, length(F))
    for u in treated_ids
        tt = t0_map[u]
        for (k, f) in enumerate(F)
            yi = df[(df.id .== u) .& (df.t .== tt + f), :y]
            y0i = df[(df.id .== u) .& (df.t .== tt + f), :y0]
            if !isempty(yi) && !isempty(y0i)
                att[k] += (yi[1] - y0i[1])
                counts[k] += 1
            end
        end
    end
    counts[counts .== 0] .= 1
    return att ./ counts
end
