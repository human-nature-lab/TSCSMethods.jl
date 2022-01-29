# balancing.jl

function balance!(model::VeryAbstractCICModel, dat)

  meanbalance!(model, dat);
  grandbalance!(model);

  return model
end
