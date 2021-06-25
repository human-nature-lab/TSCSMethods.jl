module tscsmethods

  include("dependencies.jl")

  include("pkg_types.jl")
  include("matching_functions.jl")
  include("balancing_functions.jl")
  include("estimation_functions.jl")
  include("get_wits.jl")
  include("bootstrap_estimation.jl")
  include("plotting_functions.jl")
  include("wrappers.jl")
  include("internal_helpers.jl")
  include("weekly_att.jl")

  export cicmodel,
    matching!,
    getbalance!, handle_balance!,
    estimate!,
    handle_att!, handle_attsum,
    caliper!,
    namemodel, StratDict

end

# module MyModule
#          x() = "x"
#          export x
#        end
# using .MyModule
