# getmatches!.jl

function eligiblematches!(
  observations, matches,
  rg, trtg, ids, fmin, fmax, treatcat; exg = nothing
)
  
  for (ob, tob) in zip(observations, matches)
    (tt, tu) = ob;
    @unpack mus = tob;
    if isnothing(exg)
      matchespossible!(mus, rg, trtg, ids, tt, tu, fmin, fmax, treatcat);
    else
      matchespossible!(mus, rg, trtg, ids, tt, tu, fmin, fmax, treatcat, exg);
    end
  end
  return matches
end

function matchespossible!(
  mus, rg, trtg, uid, tt, tu, fmin, fmax, treatcat::Function
)

  for (mu, rmus) in zip(uid, eachrow(mus))

    # clearly, if same unit, match is not allowed
    if tu == mu
      rmus .= false
    else

      # requires that panels are balanced
      # will need to handle missingness case eventually
      # data are padded

      # relevant treatment history for tu (treated unit)
      tu_trtimes = rg[(tt, tu)][trtg[(tt, tu)]];
      # relevant treatment history for mu (pot. match unit)
      mu_trtimes = rg[(tt, mu)][trtg[(tt, mu)]];

      fpossible!(
        rmus, tu_trtimes, mu_trtimes,
        fmin, fmax, tt,
        treatcat
      );

      # old full sliding window
      # sliding window in which matched units are barred from being treated
      # the window is [t + f - fmax, t + f - fmin]
      # if length(mtts) > 0
      #   for ttmu in mtts
      #     if tt == ttmu
      #       for α in 1:length(rmus); rmus[α] = false end
      #     elseif tt - ttmu > 0
      #       # remove f and prior fs
      #       priorremove = fmax - (tt - ttmu) # bar these
      #       for α in 1:(priorremove - fmin + 1)
      #         rmus[α] = false
      #       end
      #     elseif ttmu - tt > 0
      #       # remove f and later fs
      #       postremove = (ttmu - tt) + fmin;
      #       for α in (postremove - fmin + 1):length(rmus); rmus[α] = false end
      #     end
      #   end
      # end
      
    end
  end
  return mus
end

function matchespossible!(
  mus, rg, trtg, uid, tt, tu, fmin, fmax, treatcat::Function, exg
)

  for (mu, rmus) in zip(uid, eachrow(mus))

    # clearly, if same unit, match is not allowed
    if tu == mu
      rmus .= false
    else

      # requires that panels are balanced
      # will need to handle missingness case eventually
      # data are padded

      # relevant treatment history for tu (treated unit)
      tu_trtimes = rg[(tt, tu)][trtg[(tt, tu)]];
      # relevant treatment history for mu (pot. match unit)
      mu_trtimes = rg[(tt, mu)][trtg[(tt, mu)]];

      tu_exposures = exg[(tt, tu)][trtg[(tt, tu)]];
      mu_exposures = exg[(tt, mu)][trtg[(tt, mu)]];

      fpossible!(
        rmus, tu_trtimes, mu_trtimes,
        fmin, fmax, tt,
        treatcat,
        tu_exposures, mu_exposures
      );


      # old full sliding window
      # sliding window in which matched units are barred from being treated
      # the window is [t + f - fmax, t + f - fmin]
      # if length(mtts) > 0
      #   for ttmu in mtts
      #     if tt == ttmu
      #       for α in 1:length(rmus); rmus[α] = false end
      #     elseif tt - ttmu > 0
      #       # remove f and prior fs
      #       priorremove = fmax - (tt - ttmu) # bar these
      #       for α in 1:(priorremove - fmin + 1)
      #         rmus[α] = false
      #       end
      #     elseif ttmu - tt > 0
      #       # remove f and later fs
      #       postremove = (ttmu - tt) + fmin;
      #       for α in (postremove - fmin + 1):length(rmus); rmus[α] = false end
      #     end
      #   end
      # end
      
    end
  end
  return mus
end

"""
check possibility of the match for each f
"""
function fpossible!(
  rmus, tu_trtimes, mu_trtimes,
  fmin, fmax, tt,
  treatcat
)
  for φ in eachindex(rmus)
    f = φ + fmin - 1;
    tu_treatments = 0; # count for tu
    mu_treatments = 0; # count for mu
    
    # given that f, check possibility of match
    # by examining crossover period
    for h in 1:max(length(mu_trtimes), length(tu_trtimes))
      
      # post-treatment crossover window
      tx_l = tt; tx_u = tt + f - fmin;
      # pre-treatment crossover window
      ptx_l = (tt + f - fmax); ptx_u = (tt - 1);
      
      # bar on post-treatment match window
      if (length(mu_trtimes) > 0) & (length(tx_l:tx_u) > 0)
        for mtt in mu_trtimes
          if (mtt >= tx_l) & (mtt <= tx_u)
            rmus[φ] = false
            break
          end
        end
      end

      if !rmus[φ]
        # skip to next f if it is cancelled
        break
      end

      # (if the f wasn't blocked by post-treatment requirements)
      # match on pre-treatment xover: t+F-Fmax to t-1
      if (h <= length(tu_trtimes)) & (h > 0)
        if (tu_trtimes[h] >= ptx_l) & (tu_trtimes[h] <= ptx_u);
          tu_treatments += 1
        end
      end
      
      # match on pre-treatment xover: t+F-Fmax to t-1
      if (h <= length(mu_trtimes)) & (h > 0)
        if (mu_trtimes[h] >= ptx_l) & (mu_trtimes[h] <= ptx_u)
          mu_treatments += 1
        end
      end

      # check similarity of tu_treatments & mu_treatments
      # based on treatcat function (default or user defined)
      if treatcat(tu_treatments) != treatcat(mu_treatments)
        rmus[φ] = false
      end
    end
  end
  return rmus
end

"""
check possibility of the match for each f

version with exposure
"""
function fpossible!(
  rmus, tu_trtimes, mu_trtimes,
  fmin, fmax, tt,
  treatcat,
  tu_exposures, mu_exposures
)
  for φ in eachindex(rmus)
    f = φ + fmin - 1;
    
    # this needs to be a dictionary now
    # exposure => # times
    tu_treatments = Dict{Int, Int}(); # counts for tu
    mu_treatments = Dict{Int, Int}(); # counts for mu

    exposures = union(tu_exposures, mu_exposures);
    for e in exposures
      tu_treatments[e] = 0
      mu_treatments[e] = 0
    end

    # given that f, check possibility of match
    # by examining crossover period
    for h in 1:max(length(mu_trtimes), length(tu_trtimes))
      
      # post-treatment crossover window
      tx_l = tt; tx_u = tt + f - fmin;
      # pre-treatment crossover window
      ptx_l = (tt + f - fmax); ptx_u = (tt - 1);
      
      # bar on post-treatment match window
      if (length(mu_trtimes) > 0) & (length(tx_l:tx_u) > 0)
        for mtt in mu_trtimes
          if (mtt >= tx_l) & (mtt <= tx_u)
            rmus[φ] = false
            break
          end
        end
      end

      if !rmus[φ]
        # skip to next f if it is cancelled
        break
      end

      # (if the f wasn't blocked by post-treatment requirements)
      # match on pre-treatment xover: t+F-Fmax to t-1
      if (h <= length(tu_trtimes)) & (h > 0)
        if (tu_trtimes[h] >= ptx_l) & (tu_trtimes[h] <= ptx_u);
          tu_treatments[tu_exposures[h]] += 1
        end
      end
      
      # match on pre-treatment xover: t+F-Fmax to t-1
      if (h <= length(mu_trtimes)) & (h > 0)
        if (mu_trtimes[h] >= ptx_l) & (mu_trtimes[h] <= ptx_u)
          mu_treatments[mu_exposures[h]] += 1
        end
      end

      # check similarity of tu_treatments & mu_treatments
      # based on treatcat function (default or user defined)
      for e in exposures
        if treatcat(tu_treatments[e]) != treatcat(mu_treatments[e])
          # if the categories differ for any exposure, not eligible
          rmus[φ] = false
        end
      end
    end
  end
  return rmus
end
