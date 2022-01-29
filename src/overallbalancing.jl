# overallbalancing.jl

"""
    grandbalance!(model::AbstractCICModel)

Calculate the overall mean covariate balance for a model.
"""
function grandbalance!(model::AbstractCICModel)

  @unpack F, L, covariates, timevary = model;
  @unpack grandbalances, meanbalances = model;

  __grandbalance!(
    grandbalances, meanbalances,
    covariates, timevary, length(L)
  );

  return model
end

# not stratified
function __grandbalance!(
  grandbalances, meanbalances, covariates, timevary, Len
)
  for covar in covariates
    if timevary[covar]
      covec = meanbalances[!, covar];
      gbc = disallowmissing(_grandbalance(covec, Len))
      grandbalances[covar] = gbc
    else
      covec = meanbalances[!, covar]
      grandbalances[covar] = _grandbalance(covec)
    end
  end
  return grandbalances
end

"""
    grandbalance!(model::AbstractCICModelStratified)

Calculate the overall mean covariate balance for a model.
"""
function grandbalance!(stratmodel::AbstractCICModelStratified)

  @unpack F, L, covariates, timevary, strata = stratmodel;
  @unpack grandbalances, meanbalances = stratmodel;

  Len = length(L);

  if isempty(strata)
    println("N.B. stratum not defined in meanbalances. so, calculating nonstratified grandbalances.")
    @reset grandbalances = GrandDictNoStrat();
    __grandbalance!(
      grandbalances, meanbalances,
      covariates, timevary, length(L)
    );
  else
    us = sort(unique(strata));
    [grandbalances[s] = GrandDictNoStrat() for s in us];
    
    for s in us
      mb_s = @view meanbalances[strata .== s, :];
      __grandbalance!(
        grandbalances, mb_s,
        covariates, timevary, Len,
        s
      )
    end
  end
  return stratmodel
end

# stratified
function __grandbalance!(
  grandbalances, meanbalances, covariates, timevary, Len, s
)
  for covar in covariates
    if timevary[covar]
      gbc = disallowmissing(
        _grandbalance(
          meanbalances[!, covar],
          Len
        )
      )
      grandbalances[s][covar] = gbc
    else
      grandbalances[s][covar] = _grandbalance(
        meanbalances[!, covar]
      )
    end
  end
  return grandbalances
end

# common inner functions

function _grandbalance(covec, Len)
  reduced = reduce(vcat, covec);
  means = Vector{Union{Missing, Float64}}(undef, Len);
  lnct = 0
  for l in eachindex(1:Len)
    lnct += 1
    means[l] = mean(skipmissing([vec[l] for vec in reduced]))
  end
  return means
end

function _grandbalance(covec)
  return mean(skipmissing(reduce(vcat, covec)))
end
