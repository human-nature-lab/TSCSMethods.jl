using Parameters

#=
model object
* the data is kept separately
=#
ReforDataFrame = Union{Base.RefValue{DataFrame}, DataFrame}; # type for model
StratDict = Dict{Int64, Dict{Symbol, Float64}};

@with_kw mutable struct cicmodel
  title::String = "standard"
  id::Symbol
  t::Symbol
  outcome::Symbol
  treatment::Symbol
  matchingcovar::Vector{Symbol}
  caliper::Union{Dict{Symbol, Float64}, StratDict} = Dict{Symbol, Float64}()
  fmin::Int64
  fmax::Int64
  stratvar::Symbol = Symbol()
  tpoint::Int64
  tmin::Int64 # gives length of matching period with tpoint
  refinementnum::Int64 = Int64(0)
  
  # store either dataframe or reference to dataframe, dep on use cases
  matches::ReforDataFrame = Base.RefValue{DataFrame}()
  matches5::ReforDataFrame = Base.RefValue{DataFrame}()

  data::ReforDataFrame = Base.RefValue{DataFrame}()
  
  # balances pre refinement
  balances_pre::DataFrame = DataFrame()
  meanbalances_pre::DataFrame = DataFrame()
  # post refinement
  balances_post::DataFrame = DataFrame()
  meanbalances_post::DataFrame = DataFrame()

  # we may also want to store the full bootstrap estimates

  # one results per model, so store the DF
  boot_iterations::Int64 = 500
  results::DataFrame = DataFrame()
  boot_estimates::Array{Float64} = Array{Float64}(undef)

  # plot slots

  # covariate balance pre and post refinement
  pl_cb_pre::Plot = Plot()
  pl_cb_post::Plot = Plot()

  # att estimates pre and post refinement
  pl_att_pre::Plot = Plot()
  pl_att_post::Plot = Plot()
  
  # total number of treated observations
  treatednum::Union{Int64, Dict{Int64, Int64}} = Int64(0)
  
  # number of treated left over after filtering or caliper
  treatedleft::Union{Int64, Dict{Int64, Int64}} = Int64(0)
end
