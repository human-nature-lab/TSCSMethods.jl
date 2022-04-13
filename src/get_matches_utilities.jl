# utilities for fpossible!()

function define_xover_windows(tt, f, fmin, fmax)
    # post-treatment crossover window
    tx = (tt, tt + f - fmin)
    # pre-treatment crossover window
    ptx = (tt + f - fmax), (tt - 1);
    return tx, ptx
end

function block_postxover!(rmus, tx, mu_trtimes, φ)
    if (length(mu_trtimes) > 0) & (length(tx[1]:tx[2]) > 0)
        for mtt in mu_trtimes
            if (mtt >= tx[1]) & (mtt <= tx[2])
                rmus[φ] = false
                break
            end
        end
    end
    return rmus
end
   
function match_prexover!(
    tu_treatments, mu_treatments,
    tu_trtimes, mu_trtimes,
    ptx,
    h
)
    # (if the f wasn't blocked by post-treatment requirements)
    # match on pre-treatment xover: t+F-Fmax to t-1
    if (h <= length(tu_trtimes)) & (h > 0)
        if (tu_trtimes[h] >= ptx[1]) & (tu_trtimes[h] <= ptx[2]);
            tu_treatments += 1
        end
    end
        
    # match on pre-treatment xover: t+F-Fmax to t-1
    if (h <= length(mu_trtimes)) & (h > 0)
        if (mu_trtimes[h] >= ptx[1]) & (mu_trtimes[h] <= ptx[2])
            mu_treatments += 1
        end
    end
    return tu_treatments, mu_treatments
end

function match_prexover!(
    tu_treatments, mu_treatments,
    tu_trtimes, mu_trtimes,
    tu_exposures, mu_exposures,
    ptx,
    h
)
    # (if the f wasn't blocked by post-treatment requirements)
    # match on pre-treatment xover: t+F-Fmax to t-1
    if (h <= length(tu_trtimes)) & (h > 0)
        if (tu_trtimes[h] >= ptx[1]) & (tu_trtimes[h] <= ptx[2]);
            tu_treatments[tu_exposures[h]] += 1
        end
    end
        
    # match on pre-treatment xover: t+F-Fmax to t-1
    if (h <= length(mu_trtimes)) & (h > 0)
        if (mu_trtimes[h] >= ptx[1]) & (mu_trtimes[h] <= ptx[2])
            mu_treatments[mu_exposures[h]] += 1
        end
    end

    return tu_treatments, mu_treatments
end

function exposure_assign!(
    rmus, exposures, tu_treatments, mu_treatments, treatcat
)
    # check similarity of tu_treatments & mu_treatments
    # based on treatcat function (default or user defined)
    for e in exposures
        if treatcat(tu_treatments[e]) != treatcat(mu_treatments[e])
            # if the categories differ for any exposure, not eligible
            rmus[φ] = false
        end
    end
    return rmus
end

"""
    trim_model(model)
  
Remove treated observations that do not have any valid matches. This copies!
"""
function trim_model(model2)
  # remove treated observations with no valid mus
  anymatches = fill(true, length(model2.observations));
  for (i, e) in enumerate(model2.matches)
    anymatches[i] = any(e.mus)
  end

  model2 = @set model2.observations = model2.observations[anymatches];
  model2 = @set model2.matches = model2.matches[anymatches];
  model2 = @set model2.treatednum = length(model2.observations);
  return model2
end
