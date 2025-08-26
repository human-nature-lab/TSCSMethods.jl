# retrieve_matches!.jl

function eligiblematches!(
  observations, matches,
  rg, trtg, ids, fmin, fmax, treatcat; exg = nothing
)
  
  for (ob, tob) in zip(observations, matches)
    (tt, treated_unit) = ob;
    (; eligible_matches) = tob;

    matchespossible!(
        eligible_matches, rg, trtg, ids, tt, treated_unit, fmin, fmax,
        treatcat; exg = exg
      );
  end
  return matches
end

function matchespossible!(
  eligible_matches_matrix, rg, trtg, uid, tt, treated_unit, fmin, fmax, treatcat::Function;
  exg = nothing
)

  for (matched_unit, match_eligibility_row) in zip(uid, eachrow(eligible_matches_matrix))

    # clearly, if same unit, match is not allowed
    if treated_unit == matched_unit
      match_eligibility_row .= false
    else

      # N.B., data are padded

      # relevant treatment history for treated_unit & matched_unit
      treated_unit_trtimes = rg[(tt, treated_unit)][trtg[(tt, treated_unit)]];
      matched_unit_trtimes = rg[(tt, matched_unit)][trtg[(tt, matched_unit)]];

      if isnothing(exg)

        fpossible!(
          match_eligibility_row, treated_unit_trtimes, matched_unit_trtimes,
          fmin, fmax, tt,
          treatcat
        );
      else
        treated_unit_exposures = exg[(tt, treated_unit)][trtg[(tt, treated_unit)]];
        matched_unit_exposures = exg[(tt, matched_unit)][trtg[(tt, matched_unit)]];

        fpossible!(
          match_eligibility_row, treated_unit_trtimes, matched_unit_trtimes,
          fmin, fmax, tt,
          treatcat,
          treated_unit_exposures, matched_unit_exposures
        );
      end
        # Prior version with different computation: full sliding window
        # This was a sliding window in which matched units were barred from being treated
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
  return eligible_matches_matrix
end

"""
check possibility of the match for each f
"""
function fpossible!(
  match_eligibility_row, treated_unit_trtimes, matched_unit_trtimes,
  fmin, fmax, tt,
  treatcat
)

  for outcome_period_index in eachindex(match_eligibility_row)
    f = outcome_period_index + fmin - 1;
    treated_unit_treatments = 0; # count for treated_unit
    matched_unit_treatments = 0; # count for matched_unit
    
    # given that f, check possibility of match
    # by examining crossover period
    for h in 1:max(length(matched_unit_trtimes), length(treated_unit_trtimes))
      
      tx, ptx = define_xover_windows(tt, f, fmin, fmax)
        
      block_postxover!(match_eligibility_row, tx, matched_unit_trtimes, outcome_period_index)

      if !match_eligibility_row[outcome_period_index]
        # skip to next f if it is cancelled
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
        match_eligibility_row[outcome_period_index] = false
      end
    end
  end
  return match_eligibility_row
end

"""
check possibility of the match for each f

version with exposure
"""
function fpossible!(
  match_eligibility_row, treated_unit_trtimes, matched_unit_trtimes,
  fmin, fmax, tt,
  treatcat,
  treated_unit_exposures, matched_unit_exposures
)
  # Pre-compute exposure set once outside the loop
  exposures = union(treated_unit_exposures, matched_unit_exposures)
  
  # Pre-allocate dictionaries once with known exposures
  treated_unit_treatments = Dict{Int, Int}()
  matched_unit_treatments = Dict{Int, Int}()
  sizehint!(treated_unit_treatments, length(exposures))
  sizehint!(matched_unit_treatments, length(exposures))
  
  for outcome_period_index in eachindex(match_eligibility_row)
    f = outcome_period_index + fmin - 1;
    
    # Reset counts efficiently (reuse existing dictionaries)
    empty!(treated_unit_treatments)
    empty!(matched_unit_treatments)
    for e in exposures
      treated_unit_treatments[e] = 0
      matched_unit_treatments[e] = 0
    end

    # given that f, check possibility of match
    # by examining crossover period
    for h in 1:max(length(matched_unit_trtimes), length(treated_unit_trtimes))
      
      tx, ptx = define_xover_windows(tt, f, fmin, fmax)
        
      # bar on post-treatment match window
      block_postxover!(match_eligibility_row, tx, matched_unit_trtimes, outcome_period_index)

      if !match_eligibility_row[outcome_period_index]
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
  return match_eligibility_row
end
