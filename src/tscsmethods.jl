module tscsmethods

  include("dependencies.jl")

  include("types.jl")
  include("construction.jl")
  include("matching.jl")
  include("caliper.jl")
  include("balancing.jl")
  include("estimation.jl")
  include("stratification.jl")
  include("refine.jl")
  include("autobalancing.jl")
  include("model.jl")
  include("plotting.jl")

  export
    # types
    AbstractCICModel, AbstractCICModelStratified,
    CIC, CICStratified, CaliperCIC, CaliperCICStratified,
    RefinedCIC, RefinedCaliperCIC, RefinedCaliperCICStratified,
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
    modelrecord
end
