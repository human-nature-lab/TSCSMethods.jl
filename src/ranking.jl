# ranking.jl

@unpack matches = model;

function rank!(matches, flen)
  Threads.@threads for i in 1:length(matches)
    tob = @views matches[i];
    rankmatches!(tob, flen)
  end
  return matches
end

function rankmatches!(tob, flen)
  @unpack mus, fs, mudistances, ranks = tob;

  mamat = fill(Inf, flen, length(positions)); # in id order
  _mahaposition!(mamat, fs, mudistances, 1:flen);

  fbools = sum(fs[mus]) .> 0; # valid fs for treated obs
  _rankmatches!(ranks, fbools, mamat);
  
  return tob
end

function _mahaposition!(mamat, fs, mudistances, Φ)
  
  mcnt = 0
  for (m, mu) in enumerate(mus) # matchunit by matchunit
    if mu
      mcnt += 1 # index mudistances

      fsm = @views fs[m];

      φcnt = 0
      for (φ, fb) in enumerate(fsm)
        if fb
          φcnt += 1
          mamat[φ, m] = mudistances[mcnt][φcnt][1] # maha distance at 1
        end
      end
    end
  end
  return mamat
end

function _rankmatches!(ranks, fbools, mamat)
  for (φ, (fbool, ec)) in enumerate(zip(fbools, eachrow(mamat)))
    if fbool
      # in-place ranking
      # ranks[φ] = StatsBase.ordinalrank(ec);
      # ranks[φ][isinf.(ec)] .= 0; # make zero if not valid

      # really, we want to be able to shut off the mus easily
      # and make it such that we get the mu elements in ranked order
      # by simple indexing:
      ranks[φ] = @views mus[sortperm(ec)];

      # remove impossible mus (for that f)
      ranks[φ][findfirst(isinf.(ec[sortperm(ec)])):end] .= 0
        # Inf values should never be used, in any model
      # example refinement application:
      # ormus[1+5:end] .= 0
    end
  end
  return ranks
end
