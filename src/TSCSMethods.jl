module TSCSMethods

  funlist = [
    "dependencies.jl",
    "types.jl",
    "construction.jl",
    "matching_setup.jl",
    "groupindices.jl",
    "get_matches_utilities.jl",
    "getmatches!.jl",
    "getmatches!_missing.jl",
    "distancing_utilities.jl",
    "distancing.jl",
    "match!.jl",
    "ranking.jl",
    "caliper.jl",
    "meanbalancing.jl",
    "balancing_utilities.jl",
    "overallbalancing.jl",
    "balancing.jl",
    "estimation_setup.jl",
    "estimation_utilities.jl",
    "estimation_observationweights.jl",
    "estimation.jl",
    "estimation_stratified.jl",
    "overall.jl",
    "bayesfactor.jl",
    "resampling.jl",
    "bootstrapping.jl",
    "stratification.jl",
    "refine.jl",
    "autobalancing.jl",
    "model_utilities.jl",
    "storage.jl",
    "information.jl",
    "imputation.jl",
    "inspection.jl",
    "filterunits!.jl"
  ];

  for file in funlist
    include(file)
  end

  # vignette
  """
      example_data(; n_units = 100, n_days = 90, seed = 123)
  
  Generate synthetic example data for testing and demonstration purposes.
  
  # Arguments
  - `n_units::Int`: Number of units (counties/regions) to generate (default: 100)
  - `n_days::Int`: Number of time periods to generate (default: 90)
  - `seed::Int`: Random seed for reproducibility (default: 123)
  
  # Returns
  - `DataFrame`: Synthetic panel data with columns: date, fips, pop_dens, cumul_death_rate, death_rte, gub, day
  
  # Description
  Creates synthetic panel data mimicking the structure used in Feltham et al. (2023):
  - Units represent counties with FIPS codes
  - Treatment (gubernatorial elections) affects first 20% of units
  - Realistic population density and outcome variables with trends
  - Includes both time-invariant (pop_dens) and time-varying (cumul_death_rate) covariates
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

  export
    # types
    VeryAbstractCICModel, AbstractCICModel, AbstractCICModelStratified,
    CIC, CICStratified, CaliperCIC, CaliperCICStratified,
    RefinedCIC, RefinedCICStratified, RefinedCaliperCIC, RefinedCaliperCICStratified,
    # mechanics
    match!,
    balance!, checkbalances, autobalance,
    estimate!,
    stratify, variablestrat, combostrat, customstrat,
    makemodel,
    caliper, refine,
    # saving
    name_model,
    makerecords,
    default_treatmentcategories,
    showmatches, matchinfo, obsinfo,
    save,
    # utilities
    matchprocess, quick_att, variable_filter, treatedinfo,
    relabel!, trim_model,
    # inspection
    inspection, pretreatment,
    # vignette
    example_data, policy_data, economic_data,
    filterunits!
end
