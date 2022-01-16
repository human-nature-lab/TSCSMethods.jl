module tscsmethods

  include("dependencies.jl")

  include("types.jl")
  include("construction.jl")
  include("matching_setup.jl")
  include("makegroup_indices.jl")
  include("getmatches!.jl")
  include("distancing.jl")
  include("match!.jl")
  include("ranking.jl")
  include("caliper.jl")
  include("balancing.jl")
  include("estimation_observationweights.jl")
  include("estimation.jl")
  include("stratification.jl")
  include("refine.jl")
  include("autobalancing.jl")
  include("model.jl")
  include("plotting.jl")
  include("information.jl")

  export
    # types
    AbstractCICModel, AbstractCICModelStratified,
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
    modelrecord, save_record, save_records,
    default_treatmentcategories,
    showmatches, matchinfo, obsinfo
end
