# getmatches!.jl

function setup_missing(fmin, fmax)
  validfs = fill(true, length(fmin:fmax))
  return validfs
end

function checkmissing!(validfs, ob, rnge, rg, tg)
  for (r, γ) in zip(rg[ob], tg[ob][:, 1])

    if ismissing(γ)
      for (k, τ) in enumerate(rnge)
        if (k == 1) & (τ == r)
          # mus[:] .= false # cancel all matches
          break
          for i in eachindex(validfs)
            validfs[i] = false
          end
          return validfs
        elseif τ == r
          validfs[k-1] = false
          # mus[:, k-1] .= false # cancel matches
          break # skip to next (r, γ)
        end
      end
    end
  end
  return validfs
end

function eligiblematches!(
  observations, matches,
  tg, rg, trtg, ids, fmin, fmax, treatcat; exg = nothing
)
  
  # (ob, tob) = collect(zip(observations, matches))[1]
  for (ob, tob) in zip(observations, matches)
    (tt, tu) = ob;
    @unpack mus = tob;

    validfs = setup_missing(fmin, fmax)
    rnge = vcat(tt-1, (tt+fmin):(tt+fmax))
    checkmissing!(validfs, ob, rnge, rg, tg)
    for (i, e) in enumerate(validfs)
      if !e
        mus[:, i] .= false
      end
    end

    if any(validfs) # if the obs is possible
      matchespossible!(
        mus, rg, trtg, ids, tt, tu, fmin, fmax,
        tg, validfs, rnge,
        treatcat; exg = exg
      );
    end
  end
  return matches
end

function matchespossible!(
  mus, rg, trtg, uid, tt, tu, fmin, fmax,
  tg, validfs, rnge,
  treatcat::Function;
  exg = nothing
)

  # (mu, rmus) = collect(zip(uid, eachrow(mus)))[3]
  for (mu, rmus) in zip(uid, eachrow(mus)) # over pot. matches

    # clearly, if same unit, match is not allowed
    if tu == mu
      rmus .= false
    else

      # check missingness
      mu_validfs = setup_missing(fmin, fmax)
      checkmissing!(mu_validfs, (tt, mu), rnge, rg, tg)

      # requires that panels are balanced
      # will need to handle missingness case eventually
      # data are padded

      # relevant treatment history for tu (treated unit) & mu
      tu_trtimes = rg[(tt, tu)][trtg[(tt, tu)]];
      mu_trtimes = rg[(tt, mu)][trtg[(tt, mu)]];

      if isnothing(exg)
        fpossible_mis!(
          rmus, tu_trtimes, mu_trtimes,
          fmin, fmax, tt,
          validfs, mu_validfs,
          treatcat
        );
      else
        tu_exposures = exg[(tt, tu)][trtg[(tt, tu)]];
        mu_exposures = exg[(tt, mu)][trtg[(tt, mu)]];

        fpossible_mis!(
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
function fpossible_mis!(
  rmus, tu_trtimes, mu_trtimes,
  fmin, fmax, tt,
  validfs, mu_validfs,
  treatcat
)

  for φ in eachindex(rmus)
    # we now additionally check if that f is not valid
    # from the missingness check
    # if either the tu or mu is not valid at that f, skip
    if (mu_validfs[φ] == false) | (validfs[φ] == false)
      rmus[φ] = false
      continue
    else
    
      f = φ + fmin - 1;
      tu_treatments = 0; # count for tu
      mu_treatments = 0; # count for mu
      
      # given that f, check possibility of match
      # by examining crossover period
      for h in 1:max(length(mu_trtimes), length(tu_trtimes))
        
        tx, ptx = define_xover_windows(tt, f, fmin, fmax)
        
        block_postxover!(rmus, tx, mu_trtimes, φ)
  
        if !rmus[φ]
          # skip to next f if it is canceled
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
  end
  return rmus
end

"""
check possibility of the match for each f

version with exposure
"""
function fpossible_mis!(
  rmus, tu_trtimes, mu_trtimes,
  fmin, fmax, tt,
  validfs, mu_validfs,
  treatcat,
  tu_exposures, mu_exposures
)
  for φ in eachindex(rmus)

    # we now additionally check if that f is not valid
    # from the missingness check
    # if either the tu or mu is not valid at that f, skip
    if (mu_validfs[φ] == false) | (validfs[φ] == false)
      rmus[φ] = false
      continue
    else

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
        
        tx, ptx = define_xover_windows(tt, f, fmin, fmax)
        
        # bar on post-treatment match window
        block_postxover!(rmus, tx, mutrtimes, φ)
  
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
  end
  return rmus
end
