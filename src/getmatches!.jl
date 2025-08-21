# getmatches!.jl

function eligiblematches!(
  observations, matches,
  rg, trtg, ids, fmin, fmax, treatcat; exg = nothing
)
  
  for (ob, tob) in zip(observations, matches)
    (tt, tu) = ob;
    @unpack mus = tob;

    matchespossible!(
        mus, rg, trtg, ids, tt, tu, fmin, fmax,
        treatcat; exg = exg
      );
  end
  return matches
end

function matchespossible!(
  mus, rg, trtg, uid, tt, tu, fmin, fmax, treatcat::Function;
  exg = nothing
)

  for (mu, rmus) in zip(uid, eachrow(mus))

    # clearly, if same unit, match is not allowed
    if tu == mu
      rmus .= false
    else

      # data are padded

      # relevant treatment history for tu (treated) & mu (match)
      tu_trtimes = rg[(tt, tu)][trtg[(tt, tu)]];
      mu_trtimes = rg[(tt, mu)][trtg[(tt, mu)]];

      if isnothing(exg)

        fpossible!(
          rmus, tu_trtimes, mu_trtimes,
          fmin, fmax, tt,
          treatcat
        );
      else
        tu_exposures = exg[(tt, tu)][trtg[(tt, tu)]];
        mu_exposures = exg[(tt, mu)][trtg[(tt, mu)]];

        fpossible!(
          rmus, tu_trtimes, mu_trtimes,
          fmin, fmax, tt,
          treatcat,
          tu_exposures, mu_exposures
        );
      end
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
      
      tx, ptx = define_xover_windows(tt, f, fmin, fmax)
        
      block_postxover!(rmus, tx, mu_trtimes, φ)

      if !rmus[φ]
        # skip to next f if it is cancelled
        break
      end

      match_prexover!(
        tu_treatments, mu_treatments,
        tu_trtimes, mu_trtimes,
        ptx,
        h
      )

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
  # Pre-compute exposure set once outside the loop
  exposures = union(tu_exposures, mu_exposures)
  
  # Pre-allocate dictionaries once with known exposures
  tu_treatments = Dict{Int, Int}()
  mu_treatments = Dict{Int, Int}()
  sizehint!(tu_treatments, length(exposures))
  sizehint!(mu_treatments, length(exposures))
  
  for φ in eachindex(rmus)
    f = φ + fmin - 1;
    
    # Reset counts efficiently (reuse existing dictionaries)
    empty!(tu_treatments)
    empty!(mu_treatments)
    for e in exposures
      tu_treatments[e] = 0
      mu_treatments[e] = 0
    end

    # given that f, check possibility of match
    # by examining crossover period
    for h in 1:max(length(mu_trtimes), length(tu_trtimes))
      
      tx, ptx = define_xover_windows(tt, f, fmin, fmax)
        
      # bar on post-treatment match window
      block_postxover!(rmus, tx, mu_trtimes, φ)

      if !rmus[φ]
        # skip to next f if it is cancelled
        break
      end
      
      match_prexover!(
        tu_treatments, mu_treatments,
        tu_trtimes, mu_trtimes,
        tu_exposures, mu_exposures,
        ptx,
        h
      )
  
      exposure_assign!(
        rmus, exposures, tu_treatments, mu_treatments, treatcat
      )
    end
  end
  return rmus
end
