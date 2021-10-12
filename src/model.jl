# model.jl

"""
    save_modelset(
      filename,
      cc::Union{cicmodel, Nothing};
      ccr::Union{refinedcicmodel, Nothing} = nothing,
      cal::Union{calipercicmodel, Nothing} = nothing,
      calr::Union{refinedcicmodel, Nothing} = nothing,
      labels::Union{Dict, Nothing} = nothing
    )

Save all models together with labels. Maintains order whether or not some model exists, requires that initial unrefined non-caliper cicmodel is present.
"""
function save_modelset(
  filename,
  cc::Union{cicmodel, Nothing};
  ccr::Union{refinedcicmodel, Nothing} = nothing,
  cal::Union{calipercicmodel, Nothing} = nothing,
  calr::Union{refinedcicmodel, Nothing} = nothing,
  labels::Union{Dict, Nothing} = nothing
)

  save_object(filename * ".jld2", [cc, ccr, cal, calr, labels])
end

function load_modelset(path_to_model)
  X = load_object(path_to_model);
  # cc, ccr, cal, calr, labels
  return X[1], X[2], X[3], X[4], X[5]
end 

"""
    name_model(cc::AbstractCICModel)

Generate the filename for a set of models.
"""
function name_model(cc::AbstractCICModel)
  strat = Symbol("") == cc.stratifier ? "" : string(cc.stratifier)
  return cc.title * "_" * string(cc.outcome) * "_" * strat
end
