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
  function example_data(; n_units::Int = 100, n_days::Int = 90, seed::Int = 123)
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
    example_data,
    filterunits!
end
