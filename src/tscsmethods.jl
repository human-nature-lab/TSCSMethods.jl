module tscsmethods

  include("dependencies.jl")

  include("types.jl")
  include("construction.jl")
  include("matching.jl")
  include("caliper.jl")
  include("balancing.jl")
  include("estimation.jl")
  include("stratification.jl")
  include("model.jl")
  include("plotting.jl")

  export
    # types
    AbstractCICModel, cicmodel, calipercicmodel, refinedcicmodel,
    # mechanics
    match!,
    balance!, balancecheck,
    make_refined, make_caliper,
    estimate!,
    stratify!, variablestrat!, combostrat!, customstrat!,
    inspectcaliper,
    autobalance,
    # plotting
    plot_cb, plot_cbs,
    model_pl, plot_modelset,
    # saving
    name_model,
    save_modelset,
    load_modelset
end
