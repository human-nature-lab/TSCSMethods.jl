# retrieve_matches_missing.jl

function setup_missing(fmin, fmax)
  validfs = fill(true, length(fmin:fmax))
  return validfs
end

function checkmissing!(validfs, ob, rnge, rg, tg)
  for (r, γ) in zip(rg[ob], tg[ob][:, 1])

    if ismissing(γ)
      for (k, τ) in enumerate(rnge)
        if (k == 1) & (τ == r)
          # eligible_matches[:] .= false # cancel all matches
          break
          for i in eachindex(validfs)
            validfs[i] = false
          end
          return validfs
        elseif τ == r
          validfs[k-1] = false
          # eligible_matches[:, k-1] .= false # cancel matches
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
  
  for (ob, tob) in zip(observations, matches)
    (tt, treated_unit) = ob;
    (; eligible_matches) = tob;

    validfs = setup_missing(fmin, fmax)
    rnge = vcat(tt-1, (tt+fmin):(tt+fmax))
    checkmissing!(validfs, ob, rnge, rg, tg)
    for (i, e) in enumerate(validfs)
      if !e
        eligible_matches[:, i] .= false
      end
    end

    if any(validfs) # if the obs is possible
      matchespossible!(
        eligible_matches, rg, trtg, ids, tt, treated_unit, fmin, fmax,
        tg, validfs, rnge,
        treatcat; exg = exg
      );
    end
  end
  return matches
end

function matchespossible!(
  eligible_matches, rg, trtg, uid, tt, treated_unit, fmin, fmax,
  tg, validfs, rnge,
  treatcat::Function;
  exg = nothing
)

  for (matched_unit, match_eligibility_row) in zip(uid, eachrow(eligible_matches)) # over pot. matches

    # clearly, if same unit, match is not allowed
    if treated_unit == matched_unit
      match_eligibility_row .= false
    else

      # check missingness
      matched_unit_validfs = setup_missing(fmin, fmax)
      checkmissing!(matched_unit_validfs, (tt, matched_unit), rnge, rg, tg)

      # requires that panels are balanced
      # will need to handle missingness case eventually
      # data are padded

      # relevant treatment history for treated_unit (treated unit) & matched_unit
      treated_unit_trtimes = rg[(tt, treated_unit)][trtg[(tt, treated_unit)]];
      matched_unit_trtimes = rg[(tt, matched_unit)][trtg[(tt, matched_unit)]];

      if isnothing(exg)
        fpossible_mis!(
          match_eligibility_row, treated_unit_trtimes, matched_unit_trtimes,
          fmin, fmax, tt,
          validfs, matched_unit_validfs,
          treatcat
        );
      else
        treated_unit_exposures = exg[(tt, treated_unit)][trtg[(tt, treated_unit)]];
        matched_unit_exposures = exg[(tt, matched_unit)][trtg[(tt, matched_unit)]];

        fpossible_mis!(
          match_eligibility_row, treated_unit_trtimes, matched_unit_trtimes,
          fmin, fmax, tt,
          treatcat,
          treated_unit_exposures, matched_unit_exposures
        );
      end

      # old full sliding window
      # sliding window in which matched units are barred from being treated
      # the window is [t + f - fmax, t + f - fmin]
      # if length(mtts) > 0
      #   for ttmatched_unit in mtts
      #     if tt == ttmatched_unit
      #       for α in 1:length(match_eligibility_row); match_eligibility_row[α] = false end
      #     elseif tt - ttmatched_unit > 0
      #       # remove f and prior fs
      #       priorremove = fmax - (tt - ttmatched_unit) # bar these
      #       for α in 1:(priorremove - fmin + 1)
      #         match_eligibility_row[α] = false
      #       end
      #     elseif ttmatched_unit - tt > 0
      #       # remove f and later fs
      #       postremove = (ttmatched_unit - tt) + fmin;
      #       for α in (postremove - fmin + 1):length(match_eligibility_row); match_eligibility_row[α] = false end
      #     end
      #   end
      # end
      
    end
  end
  return eligible_matches
end

"""
check possibility of the match for each f
"""
function fpossible_mis!(
  match_eligibility_row, treated_unit_trtimes, matched_unit_trtimes,
  fmin, fmax, tt,
  validfs, matched_unit_validfs,
  treatcat
)

  for φ in eachindex(match_eligibility_row)
    # we now additionally check if that f is not valid
    # from the missingness check
    # if either the treated_unit or matched_unit is not valid at that f, skip
    if (matched_unit_validfs[φ] == false) | (validfs[φ] == false)
      match_eligibility_row[φ] = false
      continue
    else
    
      f = φ + fmin - 1;
      treated_unit_treatments = 0; # count for treated_unit
      matched_unit_treatments = 0; # count for matched_unit
      
      # given that f, check possibility of match
      # by examining crossover period
      for h in 1:max(length(matched_unit_trtimes), length(treated_unit_trtimes))
        
        tx, ptx = define_xover_windows(tt, f, fmin, fmax)
        
        block_postxover!(match_eligibility_row, tx, matched_unit_trtimes, φ)
  
        if !match_eligibility_row[φ]
          # skip to next f if it is canceled
          break
        end
  
        match_prexover!(
          treated_unit_treatments, matched_unit_treatments,
          treated_unit_trtimes, matched_unit_trtimes,
          ptx,
          h
        )
  
        # check similarity of treated_unit_treatments & matched_unit_treatments
        # based on treatcat function (default or user defined)
        if treatcat(treated_unit_treatments) != treatcat(matched_unit_treatments)
          match_eligibility_row[φ] = false
        end

      end
    end
  end
  return match_eligibility_row
end

"""
check possibility of the match for each f

version with exposure
"""
function fpossible_mis!(
  match_eligibility_row, treated_unit_trtimes, matched_unit_trtimes,
  fmin, fmax, tt,
  validfs, matched_unit_validfs,
  treatcat,
  treated_unit_exposures, matched_unit_exposures
)
  for φ in eachindex(match_eligibility_row)

    # we now additionally check if that f is not valid
    # from the missingness check
    # if either the treated_unit or matched_unit is not valid at that f, skip
    if (matched_unit_validfs[φ] == false) | (validfs[φ] == false)
      match_eligibility_row[φ] = false
      continue
    else

      f = φ + fmin - 1;
      
      # this needs to be a dictionary now
      # exposure => # times
      treated_unit_treatments = Dict{Int, Int}(); # counts for treated_unit
      matched_unit_treatments = Dict{Int, Int}(); # counts for matched_unit
  
      exposures = union(treated_unit_exposures, matched_unit_exposures);
      for e in exposures
        treated_unit_treatments[e] = 0
        matched_unit_treatments[e] = 0
      end
  
      # given that f, check possibility of match
      # by examining crossover period
      for h in 1:max(length(matched_unit_trtimes), length(treated_unit_trtimes))
        
        tx, ptx = define_xover_windows(tt, f, fmin, fmax)
        
        # bar on post-treatment match window
        block_postxover!(match_eligibility_row, tx, matched_unit_trtimes, φ)
  
        if !match_eligibility_row[φ]
          # skip to next f if it is cancelled
          break
        end
  
        match_prexover!(
          treated_unit_treatments, matched_unit_treatments,
          treated_unit_trtimes, matched_unit_trtimes,
          treated_unit_exposures, matched_unit_exposures,
          ptx,
          h
        )
  
        exposure_assign!(
          match_eligibility_row, exposures, treated_unit_treatments, matched_unit_treatments, treatcat
        )

      end
    end
  end
  return match_eligibility_row
end
