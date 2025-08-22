# example_data.jl

# vignette data

"""
    example_data()

Load pre-existing example data for testing and demonstration.

# Returns
- `DataFrame`: Panel data with columns: date, fips, pop_dens, cumul_death_rate, death_rte, gub, day

# Description
Loads example panel data from the package's vignette directory. This is real data
from the package examples, providing a consistent dataset for testing and tutorials.

For generating synthetic data with custom parameters, use `example_data_generated()`,
`policy_data()`, or `economic_data()` instead.
"""
function example_data()
    # Load example data from CSV file
    csv_path = joinpath(@__DIR__, "..", "vignette", "example_data.csv")
    data = CSV.read(csv_path, DataFrame)
    return data
    end

    # Keep old function for backward compatibility if needed
    function example_data_generated(; n_units::Int = 100, n_days::Int = 90, seed::Int = 123)
    Random.seed!(seed)

    # Date range (similar to original data)
    start_date = Date(2021, 10, 1)
    dates = start_date:Day(1):(start_date + Day(n_days - 1))

    # Create synthetic panel data
    data = DataFrame()

    # Treatment assignment: first 20% of units are treated (gubernatorial election)
    n_treated = max(1, round(Int, n_units * 0.2))

    for unit in 1:n_units
        is_treated_unit = (unit <= n_treated)
        
        # Unit-specific characteristics
        base_pop_dens = 20.0 + 100.0 * rand()  # Population density: 20-120
        base_death_rate = 15.0 + 10.0 * rand()  # Base death rate: 15-25
        
        # Treatment effect (if treated)
        treatment_effect = is_treated_unit ? -0.15 * rand() : 0.0  # Small negative effect
        
        for (day_idx, date) in enumerate(dates)
        # Time-varying cumulative death rate with trend
        time_trend = day_idx * 0.05
        cumul_death_rate = base_death_rate + time_trend + randn() * 2.0
        cumul_death_rate = max(0.0, cumul_death_rate)  # Non-negative
        
        # Daily death rate (outcome variable)
        seasonal_effect = 0.1 * sin(2π * day_idx / 30)  # Monthly seasonality
        base_outcome = 0.5 + seasonal_effect + randn() * 0.3
        
        # Add treatment effect for treated units after day 30
        treated_outcome = (is_treated_unit && day_idx > 30) ? treatment_effect : 0.0
        death_rte = max(0.0, base_outcome + treated_outcome)
        
        # Government treatment indicator (binary) - only 1 on treatment event day
        # For staggered design: treated units get treatment at a specific random day
        treatment_day = is_treated_unit ? max(1, round(Int, n_days * 0.4)) + (unit % 10) : -1
        gub = (is_treated_unit && day_idx == treatment_day) ? 1 : 0
        
        push!(data, (
            date = date,
            fips = 1000 + unit,  # FIPS-like codes starting at 1001
            pop_dens = base_pop_dens + randn() * 2.0,  # Small time variation
            cumul_death_rate = cumul_death_rate,
            death_rte = death_rte,
            gub = gub
        ))
        end
    end

    # Add integer day variable (0-based from start date)
    data[!, :day] = repeat(0:(n_days-1), n_units)

    return data
end

# data_generation.jl
# Synthetic data generators for statistical correctness testing

using Random
using DataFrames
using Dates

"""
    generate_simple_tscs(; kwargs...)

Generate synthetic time-series cross-sectional data with known treatment effect.

# Mathematical Model
Y_it = α_i + λ_t + τ*D_it + ε_it

Where:
- Y_it: outcome for unit i at time t  
- α_i: unit fixed effects (time-invariant differences between units)
- λ_t: time fixed effects (common time trends)
- D_it: treatment indicator (1 if unit i treated at time t, 0 otherwise)
- τ: treatment effect (THIS IS WHAT WE WANT TO RECOVER)
- ε_it: idiosyncratic error term

# Key Properties
- **Parallel trends satisfied by construction** (common λ_t)
- **No confounding** (treatment assignment independent of potential outcomes)
- **Known counterfactual** (Y_it without treatment = α_i + λ_t + ε_it)

# Parameters
- `true_att::Float64`: The true average treatment effect (τ)
- `n_units::Int=100`: Number of cross-sectional units
- `n_periods::Int=50`: Number of time periods
- `treatment_period::Int=25`: When treatment begins (must be > max(abs(L)))
- `n_treated::Int=20`: Number of units that receive treatment
- `unit_fe_var::Float64=2.0`: Variance of unit fixed effects (α_i)
- `time_fe_var::Float64=1.0`: Variance of time fixed effects (λ_t)
- `error_var::Float64=1.0`: Variance of idiosyncratic errors (ε_it)
- `seed::Int=1234`: Random seed for reproducibility

# Returns
DataFrame with columns:
- `unit_id`: Unit identifier (1, 2, ..., n_units)
- `time_period`: Time identifier (1, 2, ..., n_periods)
- `treatment`: Treatment indicator (0/1)
- `outcome`: Generated outcome variable Y_it
- `unit_fe`: Unit fixed effect α_i (for validation)
- `time_fe`: Time fixed effect λ_t (for validation)
- `true_control_outcome`: Counterfactual outcome without treatment

# Example
```julia
# Generate data where true ATT = 2.5
data = generate_simple_tscs(true_att=2.5, n_units=50, n_periods=40, seed=42)

# The package should recover ATT ≈ 2.5
model = makemodel(data, :time_period, :unit_id, :treatment, :outcome, [], Dict(), 1:5, -10:-1)
match!(model, data)
estimate!(model, data)
@test abs(model.overall.ATT - 2.5) < 0.2  # Should be close to true effect
```
"""
function generate_simple_tscs(;
    true_att::Float64,
    n_units::Int = 100,
    n_periods::Int = 50,
    treatment_period::Int = 25,
    treatment_duration::Int = 10,  # NEW: Duration of treatment
    n_treated::Int = 20,
    unit_fe_var::Float64 = 2.0,
    time_fe_var::Float64 = 1.0,
    error_var::Float64 = 1.0,
    seed::Int = 1234
)
    
    # Input validation
    if treatment_period >= n_periods
        throw(ArgumentError("treatment_period ($treatment_period) must be < n_periods ($n_periods)"))
    end
    if (treatment_period + treatment_duration + 15) > n_periods  # Need at least 15 periods after treatment for F windows
        throw(ArgumentError("treatment_period ($treatment_period) + treatment_duration ($treatment_duration) + 15 must be <= n_periods ($n_periods) for adequate post-treatment periods for F windows"))
    end
    if treatment_period <= 30  # Need enough pre-treatment periods for L windows (up to -30)
        throw(ArgumentError("treatment_period ($treatment_period) must be > 30 for adequate pre-treatment periods for L windows"))
    end
    if n_treated >= n_units
        throw(ArgumentError("n_treated ($n_treated) must be < n_units ($n_units)"))
    end
    
    Random.seed!(seed)
    
    # Generate unit fixed effects α_i ~ N(0, unit_fe_var)
    unit_effects = randn(n_units) * sqrt(unit_fe_var)
    
    # Generate time fixed effects λ_t ~ N(0, time_fe_var)
    # Add slight linear trend to make more realistic
    time_effects = randn(n_periods) * sqrt(time_fe_var) + 0.02 * (1:n_periods)
    
    # Randomly assign treatment (first n_treated units get treated)
    treated_units = sort(ordered_sample(1:n_units, n_treated))
    
    # Create the dataset
    data = DataFrame()
    
    for unit in 1:n_units
        is_treated = unit in treated_units
        
        for period in 1:n_periods
            # Treatment indicator: 1 if treated unit and in treatment window
            in_treatment_window = period >= treatment_period && period < (treatment_period + treatment_duration)
            treatment_status = (is_treated && in_treatment_window) ? 1 : 0
            
            # Generate idiosyncratic error
            error_term = randn() * sqrt(error_var)
            
            # Compute outcome: Y_it = α_i + λ_t + τ*D_it + ε_it
            true_control_outcome = unit_effects[unit] + time_effects[period] + error_term
            actual_outcome = true_control_outcome + (treatment_status * true_att)
            
            push!(data, (
                unit_id = unit,
                time_period = period,
                treatment = treatment_status,
                outcome = actual_outcome,
                unit_fe = unit_effects[unit],
                time_fe = time_effects[period],
                true_control_outcome = true_control_outcome
            ))
        end
    end
    
    return data
end

"""
    generate_tscs_with_covariates(; kwargs...)

Generate synthetic TSCS data with covariates that affect both treatment and outcome.

# Mathematical Model
Y_it = α_i + λ_t + β*X_it + τ*D_it + ε_it
D_it = treatment assignment potentially related to X_it

Where X_it are time-varying covariates that may confound the relationship.

# Additional Parameters
- `covariate_effects::Vector{Float64}`: Effects of each covariate on outcome (β)
- `covariate_names::Vector{Symbol}`: Names for the covariates
- `confounding_strength::Float64=0.0`: How much covariates affect treatment assignment
"""
function generate_tscs_with_covariates(;
    true_att::Float64,
    covariate_effects::Vector{Float64} = [1.5, -0.8],
    covariate_names::Vector{Symbol} = [:covar1, :covar2],
    confounding_strength::Float64 = 0.0,
    kwargs...  # Pass through other parameters to generate_simple_tscs
)
    
    # Generate base data without confounding
    base_data = generate_simple_tscs(; true_att=true_att, kwargs...)
    
    n_units = maximum(base_data.unit_id)
    n_periods = maximum(base_data.time_period)
    n_covariates = length(covariate_effects)
    
    if length(covariate_names) != n_covariates
        throw(ArgumentError("Length of covariate_names must match covariate_effects"))
    end
    
    # Generate covariates for each observation
    Random.seed!(kwargs[:seed] + 1000)  # Different seed for covariates
    
    covariate_data = DataFrame()
    
    for i in 1:nrow(base_data)
        row_data = Dict{Symbol, Any}()
        
        # Copy existing columns
        for col in names(base_data)
            row_data[Symbol(col)] = base_data[i, col]
        end
        
        # Generate covariates (could be time-varying or unit-specific)
        unit_id = base_data[i, :unit_id]
        time_period = base_data[i, :time_period]
        
        covariate_contribution = 0.0
        for (j, (effect, name)) in enumerate(zip(covariate_effects, covariate_names))
            # Generate covariate value (mix of unit and time variation)
            covar_value = randn() + 0.3 * unit_id/n_units + 0.1 * time_period/n_periods
            row_data[name] = covar_value
            covariate_contribution += effect * covar_value
        end
        
        # Update outcome to include covariate effects
        row_data[:outcome] = base_data[i, :outcome] + covariate_contribution
        row_data[:true_control_outcome] = base_data[i, :true_control_outcome] + covariate_contribution
        
        # If confounding requested, adjust treatment assignment (complex - implement later)
        if confounding_strength > 0.0
            @warn "Confounding not yet implemented - using random assignment"
        end
        
        push!(covariate_data, row_data)
    end
    
    return covariate_data
end

"""
    validate_dgp(data::DataFrame, true_att::Float64)

Validate that the data generating process is working correctly.
This is a helper function for debugging DGPs.
"""
function validate_dgp(data::DataFrame, true_att::Float64)
    println("=== DGP Validation ===")
    
    # Check basic structure
    println("Data dimensions: $(nrow(data)) rows, $(ncol(data)) columns")
    println("Units: $(length(unique(data.unit_id)))")
    println("Time periods: $(length(unique(data.time_period)))")
    println("Treated observations: $(sum(data.treatment))")
    
    # Check that treatment effect is approximately correct
    treated_outcomes = data[data.treatment .== 1, :outcome]
    control_outcomes = data[data.treatment .== 0, :outcome]
    naive_diff = mean(treated_outcomes) - mean(control_outcomes)
    println("Naive treatment effect: $(round(naive_diff, digits=3)) (true: $true_att)")
    
    # Check parallel trends assumption (should be satisfied by construction)
    pre_treatment_data = data[data.time_period .< 25, :]  # Assuming treatment at period 25
    if nrow(pre_treatment_data) > 0
        treated_units = unique(data[data.treatment .== 1, :unit_id])
        control_units = unique(data[data.treatment .== 0, :unit_id])
        
        treated_pre = pre_treatment_data[in.(pre_treatment_data.unit_id, Ref(treated_units)), :]
        control_pre = pre_treatment_data[in.(pre_treatment_data.unit_id, Ref(control_units)), :]
        
        if nrow(treated_pre) > 0 && nrow(control_pre) > 0
            treated_pre_mean = mean(treated_pre.outcome)
            control_pre_mean = mean(control_pre.outcome)
            pre_diff = treated_pre_mean - control_pre_mean
            println("Pre-treatment difference: $(round(pre_diff, digits=3)) (should be small)")
        end
    end
    
    # Validate fixed effects structure
    if :unit_fe in names(data) && :time_fe in names(data)
        reconstructed = data.unit_fe + data.time_fe + (data.treatment * true_att)
        residuals = data.outcome - reconstructed
        println("Residual std dev: $(round(std(residuals), digits=3)) (should ≈ error_var)")
    end
    
    println("======================")
end

"""
    generate_realistic_tscs(; kwargs...)

Generate synthetic data that mimics the WORKING example_data() treatment pattern
but with known, controllable treatment effects for validation.

# Key Design
- **Event-based treatment**: Treatment indicator = 1 only on treatment day
- **Sparse treatment**: ~2% of observations treated (like real example)
- **Staggered timing**: Different units treated on different days
- **Known ATT**: Specified treatment effect we can validate against

# Parameters
- `true_att::Float64`: Known treatment effect to recover
- `n_units::Int=100`: Number of units (like example data)
- `n_days::Int=90`: Number of time periods (like example data)
- `n_treated::Int=20`: Number of treated units (like example data)
- `treatment_start_day::Int=40`: Earliest treatment day (like example data)
- `outcome_scale::Float64=0.5`: Base outcome level (like death rates)
- `seed::Int=1234`: Random seed

# Returns
DataFrame with same structure as example_data():
- `date`, `fips`, `pop_dens`, `cumul_death_rate`, `death_rte`, `gub`, `day`
- Treatment effect built into `death_rte` outcome
"""
function generate_realistic_tscs(;
    true_att::Float64,
    n_units::Int = 100,
    n_days::Int = 90,
    n_treated::Int = 20,
    treatment_start_day::Int = 40,
    outcome_scale::Float64 = 0.5,
    seed::Int = 1234
)
    Random.seed!(seed)
    
    # Mimic example_data() structure exactly
    start_date = Date(2021, 10, 1)
    dates = start_date:Day(1):(start_date + Day(n_days - 1))
    
    data = DataFrame()
    
    # Treatment assignment: first n_treated units get treated (like example)
    treated_units = 1:n_treated
    
    for unit in 1:n_units
        is_treated_unit = (unit in treated_units)
        
        # Unit-specific characteristics (like example_data)
        base_pop_dens = 20.0 + 100.0 * rand()  # Population density: 20-120
        base_death_rate = 15.0 + 10.0 * rand()  # Base death rate: 15-25
        
        # Treatment timing: staggered like example_data
        treatment_day = is_treated_unit ? treatment_start_day + (unit % 10) : -1
        
        for (day_idx, date) in enumerate(dates)
            # Time-varying cumulative death rate (like example)
            time_trend = day_idx * 0.05
            cumul_death_rate = base_death_rate + time_trend + randn() * 2.0
            cumul_death_rate = max(0.0, cumul_death_rate)
            
            # Base outcome with seasonality (like example)  
            seasonal_effect = 0.1 * sin(2π * day_idx / 30)
            base_outcome = outcome_scale + seasonal_effect + randn() * 0.3
            
            # Treatment indicator: 1 ONLY on treatment day (key pattern!)
            gub = (is_treated_unit && day_idx == treatment_day) ? 1 : 0
            
            # Treatment effect: Apply AFTER treatment day (like example)
            # This creates the ATT we want to recover
            if is_treated_unit && day_idx > treatment_day
                treated_effect = true_att
                death_rte = max(0.0, base_outcome + true_att)
            else
                treated_effect = 0.0
                death_rte = max(0.0, base_outcome)
            end
            
            push!(data, (
                date = date,
                fips = 1000 + unit,  # FIPS codes starting at 1001
                pop_dens = base_pop_dens + randn() * 2.0,  # Small time variation
                cumul_death_rate = cumul_death_rate,
                death_rte = death_rte,
                gub = gub,
                # Store the true counterfactual for validation
                true_counterfactual = base_outcome,
                true_treatment_effect = treated_effect
            ))
        end
    end
    
    # Add integer day variable (0-based like example)
    data[!, :day] = repeat(0:(n_days-1), n_units)
    
    return data
end

# Helper function for sampling without replacement
ordered_sample(x, n) = x[sort(randperm(length(x))[1:n])]

"""
    policy_data(; n_units = 30, n_days = 120, seed = 456)

Generate synthetic policy intervention dataset for demonstration.

# Arguments
- `n_units::Int`: Number of units (states/regions) to generate (default: 30)
- `n_days::Int`: Number of time periods to generate (default: 120) 
- `seed::Int`: Random seed for reproducibility (default: 456)

# Returns
- `DataFrame`: Synthetic panel data with policy intervention theme

# Description
Simulates a policy intervention study where different units adopt policies at different times.
Includes multiple time-varying covariates and realistic policy effect patterns.
"""
function policy_data(; n_units::Int = 30, n_days::Int = 120, seed::Int = 456)
    Random.seed!(seed)

    start_date = Date(2020, 1, 1)
    dates = start_date:Day(1):(start_date + Day(n_days - 1))

    data = DataFrame()

    # More units get treated (40% adoption)
    n_treated = round(Int, n_units * 0.4)

    for unit in 1:n_units
        is_treated_unit = (unit <= n_treated)
        
        # Unit characteristics  
        base_gdp = 25.0 + 15.0 * rand()  # GDP per capita (thousands)
        base_education = 0.6 + 0.3 * rand()  # Education index
        base_outcome = 50.0 + 20.0 * rand()  # Policy outcome baseline
        
        # Treatment timing - spread across middle period
        treatment_day = is_treated_unit ? round(Int, 40 + rand() * 40) : -1
        
        for (day_idx, date) in enumerate(dates)
        # Time-varying covariates with realistic trends
        gdp_trend = day_idx * 0.02 + randn() * 0.5
        gdp_percap = max(10.0, base_gdp + gdp_trend)
        
        education_trend = day_idx * 0.001 + randn() * 0.02  
        education_idx = min(1.0, max(0.0, base_education + education_trend))
        
        # Policy outcome with treatment effect
        seasonal = 2.0 * sin(2π * day_idx / 30)
        time_trend = day_idx * 0.05
        base_policy_outcome = base_outcome + seasonal + time_trend + randn() * 3.0
        
        # Treatment effect kicks in gradually after policy
        if is_treated_unit && day_idx > treatment_day
            periods_since = day_idx - treatment_day
            effect_size = 5.0 * (1 - exp(-periods_since / 10.0))  # Gradual effect
            base_policy_outcome += effect_size
        end
        
        # Policy intervention indicator  
        policy = (is_treated_unit && day_idx == treatment_day) ? 1 : 0
        
        push!(data, (
            date = date,
            state_id = 2000 + unit,
            gdp_percap = gdp_percap,
            education_idx = education_idx,  
            policy_outcome = max(0.0, base_policy_outcome),
            policy = policy,
            day = day_idx - 1
        ))
        end
    end

return data
end

"""
    economic_data(; n_units = 40, n_days = 90, seed = 789)

Generate synthetic economic shock dataset for demonstration.

# Arguments  
- `n_units::Int`: Number of units (countries/regions) to generate (default: 40)
- `n_days::Int`: Number of time periods to generate (default: 90)
- `seed::Int`: Random seed for reproducibility (default: 789)

# Returns
- `DataFrame`: Synthetic panel data with economic shock theme

# Description
Simulates an economic analysis where units experience economic shocks (e.g., natural disasters,
financial crises) at different times. Includes realistic economic indicators and shock effects.
"""
function economic_data(; n_units::Int = 40, n_days::Int = 90, seed::Int = 789)
    Random.seed!(seed)

    start_date = Date(2022, 1, 1)  
    dates = start_date:Day(7):(start_date + Day(7 * (n_days - 1)))  # Weekly data

    data = DataFrame()

    # 25% of units experience economic shock  
    n_shocked = round(Int, n_units * 0.25)

    for unit in 1:n_units
        experiences_shock = (unit <= n_shocked)
        
        # Unit economic characteristics
        base_unemployment = 0.03 + 0.07 * rand()  # 3-10% unemployment
        base_inflation = 0.01 + 0.04 * rand()     # 1-5% inflation
        base_gdp_growth = 0.005 + 0.015 * rand()  # 0.5-2% quarterly growth
        
        # Shock timing - happens in middle third of time period
        shock_day = experiences_shock ? round(Int, n_days/3 + rand() * n_days/3) : -1
        
        for (day_idx, date) in enumerate(dates)
        # Economic indicators with realistic dynamics
        cycle = 0.5 * sin(2π * day_idx / 52)  # Annual cycle
        
        unemployment_change = cycle * 0.01 + randn() * 0.005
        unemployment = max(0.01, min(0.20, base_unemployment + unemployment_change))
        
        inflation_change = cycle * 0.005 + randn() * 0.003  
        inflation = max(-0.02, min(0.15, base_inflation + inflation_change))
        
        gdp_change = randn() * 0.01
        gdp_growth = base_gdp_growth + gdp_change
        
        # Economic shock effect
        if experiences_shock && day_idx >= shock_day
            periods_since = day_idx - shock_day
            # Sharp initial impact that gradually recovers
            shock_intensity = 3.0 * exp(-periods_since / 15.0)
            
            unemployment += shock_intensity * 0.02  # Unemployment spikes
            inflation += shock_intensity * 0.01     # Inflation rises  
            gdp_growth -= shock_intensity * 0.005   # Growth drops
        end
        
        # Economic shock indicator
        econ_shock = (experiences_shock && day_idx == shock_day) ? 1 : 0
        
        push!(data, (
            date = date,
            country_id = 3000 + unit,
            unemployment = unemployment,
            inflation = inflation, 
            gdp_growth = gdp_growth,
            econ_shock = econ_shock,
            week = day_idx - 1
        ))
        end
    end

    return data
end

