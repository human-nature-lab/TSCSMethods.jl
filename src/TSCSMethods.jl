module TSCSMethods

  funlist = [
    "dependencies.jl",
    "types.jl",
    "construction.jl",
    "matching_setup.jl",
    "make_groupindices.jl",
    "get_matches_utilities.jl",
    "getmatches!.jl",
    "getmatches!_missing.jl",
    "distancing_utilitites.jl",
    "distancing.jl",
    "match!.jl",
    "ranking.jl",
    "caliper.jl",
    "meanbalancing.jl",
    "balancing_utilities.jl",
    "overallbalancing.jl",
    "balancing.jl",
    "estimation_observationweights.jl",
    "resampling.jl",
    "bootstrapping.jl",
    "estimation_setup.jl",
    "estimation.jl",
    "stratification.jl",
    "refine.jl",
    "autobalancing.jl",
    "model.jl",
    "plotting.jl",
    "information.jl",
  ];

  for file in funlist
    include(file)
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
    # plotting
    plot_cb, plot_cbs,
    model_pl, plot_modelset,
    # saving
    name_model,
    makerecords,
    default_treatmentcategories,
    showmatches, matchinfo, obsinfo,
    # utilities
    matchprocess, quick_att, variable_filter, treatedinfo,
    relabel!, trim_model
end
