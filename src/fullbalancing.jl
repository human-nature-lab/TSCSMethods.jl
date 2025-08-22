# fullbalancing.jl

"""
    fullbalance!(model::AbstractCICModel, dat::DataFrame)

calculate the full set of standardized balance scores, for each treated observation - matched unit pair.

cf. meanbalance for the averages by treated observation
"""
function fullbalance(model::AbstractCICModel, dat::DataFrame)

  (; observations, matches, ids) = model;
  (; F, L, id, t, reference, treatment, covariates, timevary) = model;
  fmin = minimum(F); fmax = maximum(F); mmin = minimum(L); mmax = maximum(L);
  
  tmin = minimum(dat[!, t])
  X = Matrix(dat[!, covariates]);

  balances = setup_fullbalances(model);

  # this function needs to be rewritten
  # balances = allocate_balances(
  #   matches, covariates, timevary, matches, ids, Lrnge
  # );
  
  Lσ = compute_treated_std(model, dat);

  tg, rg, trtg = make_groupindices(
    dat[!, t], dat[!, treatment],
    dat[!, id], ids,
    fmax, mmin,
    X
  );

  _balance!(
    balances,
    observations,
    matches,
    ids,
    fmin, fmax, mmin, mmax,
    tmin,
    tg, rg,
    covariates, timevary,
    Lσ, reference
  );

  return balances
end

function _balance!(
  balances,
  observations,
  tobs,
  ids,
  fmin, fmax, mmin, mmax,
  tmin,
  tg, rg,
  covariates, timevary,
  Lσ, reference
)

  @inbounds Threads.@threads :greedy for i in eachindex(observations)
    balrwi = @view balances[i, :];
    ob = observations[i];
    emus = matchassignments(tobs[i], ids; returnefsets = false);

    if length(emus) == 0
      continue
    end

    # we need the difference between the same covariates at each time point
    # (don't need each row)
    # we know which fs are allowable already, so we don't need to worry about lengths, etc. (since the panels are necessarily balanced)

    # fw = matchwindow(window_index + fmin - 1, tt, mmin, mmax);
    ttmatch = mmin + ob[1] + fmin : mmax + ob[1] + fmax;

    γcs = eachcol(tg[ob]);

    # treated_covariate_rows = eachrow(tg[ob]);
    # lag_times = rg[ob];
    firstγtime = rg[ob][begin];

    # if the data vector is shorter than ttmatch
    # because of late data start
    # can ignore early end case (there won't be a
    # treatment known past the data...)
    tadjustment = if minimum(ttmatch) < tmin
      abs(minimum(ttmatch))
    else 0
    end

    balance_atreated!(
      balrwi, emus, ob[1], tg, covariates, timevary, γcs, Lσ,
      tadjustment, reference, firstγtime
    );

  end
  return tobs
end

function _timevarybalance!(
  distlocs, γc, gc, Lσ, tt, covar, tadjustment, xtfirst
)
  for (l, (x, y)) in enumerate(zip(γc, gc))
    xtfirst += 1
    distlocs[l+tadjustment] = (x - y) * Lσ[xtfirst - tt, covar]
  end
end

function balance_amatch!(
  balrwi, covariates, timevary,
  γcs, gcs, Lσ, tadjustment, reference, firstγtime, tt, em
)
  for (covar, γc, gc) in zip(covariates, γcs, gcs)
    # N.B. eachrow(balances)[i][covar]
    balic = balrwi[covar]
    if timevary[covar]
      _timevarybalance!(
        balic[em], γc, gc, Lσ, tt, covar, tadjustment, firstγtime
      )
    else
      balic[em] = (γc[1] - gc[1]) * Lσ[reference, covar]
      continue
    end
  end
  return balrwi
end

function balance_atreated!(
  balrwi, emus, tt, tg, covariates, timevary, γcs, Lσ,
  tadjustment, reference, firstγtime
)
  for (em, emu) in enumerate(emus)
    # emu = emus[e];
    # efs = efsets[em]; # logical rep. of fs that exist for a match    
    gcs = eachcol(tg[(tt, emu)]);

    # do all 80, for each covar, here
    balance_amatch!(
      balrwi, covariates, timevary,
      γcs, gcs, Lσ, tadjustment, reference, firstγtime, tt, em
    )
  end
  return balrwi
end

function setup_fullbalances(model)
  (; F, L, covariates, timevary, matches, ids, observations) = model;

  Lrnge = length((F[1] + L[1]):(F[end] + L[end]));

  # row for each treated obs, with vector of treated observations
  # prob don't need observations, matchunits in dataframe?
  balances = DataFrame(
    # tob = observations
  );

  timevar = Vector{Symbol}();
  staticvar = Vector{Symbol}();
  for covar in covariates
    if timevary[covar]
      push!(timevar, covar)
      balances[!, covar] = [
        Vector{Union{Missing, Float64}}(missing, Lrnge) for _ in 1:length(observations)
      ]
      balances[!, covar] = Vector{
        Vector{Vector{Union{Missing, Float64}}}
      }(undef, length(observations));
      for i in eachindex(observations)
        balances[i, covar] = [
          Vector{Union{Missing, Float64}}(missing, Lrnge) for _ in 1:sum(matches[i].eligible_matches)
        ]
      end
    else
    push!(staticvar, covar)
    balances[!, covar] = zeros(Float64, nrow(balances));
    balances[!, covar] = Vector{Vector{Float64}}(undef, length(observations));
    for i in eachindex(observations)
      balances[i, covar] = Vector{Float64}(undef, sum(matches[i].eligible_matches))
      end
    end
  end

  return balances
end

# function matchassignments(tobsi, ids; returnefsets = true)
#   assigned = getassigned(tobsi.eligible_matches, tobsi.fs);
#   emus = ids[assigned]; # eligible matches
  
#   if returnefsets == false
#     return emus
#   else
#     efsets = tobsi.fs[assigned]; # allowable fs for each eligible match
#     return emus, efsets
#   end
# end
