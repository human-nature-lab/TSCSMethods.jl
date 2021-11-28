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
    modelrecord, save_record, save_records

  # module fullbalancing
  #   using DataFrames, DataFramesMeta
  #   using tscsmethods:AbstractCICModel, @unpack, std_treated
  #   using tscsmethods:allocate_meanbalances!
  #   using tscsmethods:make_groupindices, mean, matchassignments
  #   using tscsmethods:grandbalance!
  #   include("fullbalancing.jl")
  #   include("mean_fullbalancing.jl")

  #   export fullbalance, balances!
  # end
end
