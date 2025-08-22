# overallbalancing.jl

"""
    grandbalance!(model::AbstractCICModel)

Calculate the overall mean covariate balance for a model.
"""
function grandbalance!(model::AbstractCICModel)

  (; L, covariates, timevary, grandbalances, meanbalances) = model

  __compute_grand_balance!(
    grandbalances, meanbalances,
    covariates, timevary, length(L)
  );

  return model
end

# not stratified
function __compute_grand_balance!(
  grandbalances, meanbalances, covariates, timevary, Len
)
  for covar in covariates
    if timevary[covar]
      covariate_balance_vector = meanbalances[!, covar];
      gbc = disallowmissing(_compute_grand_balance(covariate_balance_vector, Len))
      grandbalances[covar] = gbc
    else
      covariate_balance_vector = meanbalances[!, covar]
      grandbalances[covar] = _compute_grand_balance(covariate_balance_vector)
    end
  end
  return grandbalances
end

"""
    grandbalance!(model::AbstractCICModelStratified)

Calculate the overall mean covariate balance for a model.
"""
function grandbalance!(stratmodel::AbstractCICModelStratified)

  (; L, covariates, timevary, strata, grandbalances, meanbalances) = stratmodel

  Len = length(L);

  if isempty(strata)
    println("N.B. stratum not defined in meanbalances. so, calculating nonstratified grandbalances.")
    @reset grandbalances = GrandDictNoStrat();
    __compute_grand_balance!(
      grandbalances, meanbalances,
      covariates, timevary, Len
    );
  else
    us = sort(unique(strata));
    [grandbalances[s] = GrandDictNoStrat() for s in us];
    
    for s in us
      mb_s = @view meanbalances[strata .== s, :];
      __compute_grand_balance!(
        grandbalances, mb_s,
        covariates, timevary, Len,
        s
      )
    end
  end
  return stratmodel
end

# stratified
function __compute_grand_balance!(
  grandbalances, meanbalances, covariates, timevary, Len, s
)
  for covar in covariates
    if timevary[covar]
      gbc = disallowmissing(
        _compute_grand_balance(
          meanbalances[!, covar],
          Len
        )
      )
      grandbalances[s][covar] = gbc
    else
      grandbalances[s][covar] = _compute_grand_balance(
        meanbalances[!, covar]
      )
    end
  end
  return grandbalances
end

# common inner functions

function _compute_grand_balance(covariate_balance_vector, Len)
  reduced = reduce(vcat, covariate_balance_vector);
  # means = Vector{Union{Missing, Float64}}(undef, Len);
  means = Vector{Float64}(undef, Len);

  # handling missing periods
  # where no units have a value within the matching period
  # output an NaN
  # for now, this is fine, but better to use missing,
  # and then convert to NaN where needed (e.g. Makie plotting)

  for l in eachindex(1:Len)
    vecs = [vec[l] for vec in reduced]
    if !all(ismissing.(vecs))
      means[l] = mean(skipmissing(vecs))
    else 
      means[l] = NaN
    end
  end
  return means
end

function _compute_grand_balance(covariate_balance_vector)
  # Handle BalanceData specially
  if eltype(covariate_balance_vector) <: BalanceData
    all_values = Float64[]
    for bd in covariate_balance_vector
      valid_indices = .!bd.is_missing
      append!(all_values, bd.values[valid_indices])
    end
    return isempty(all_values) ? NaN : mean(all_values)
  else
    return mean(skipmissing(reduce(vcat, covariate_balance_vector)))
  end
end
