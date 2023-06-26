# ranking.jl

## ranking for the sliding window

function rank!(matches, flen)
  Threads.@threads for i in 1:length(matches)
    tob = @views matches[i];
    rankmatches!(tob, flen)
  end
  return matches
end

function rankmatches!(tob, flen)
  @unpack mus, distances, ranks = tob;

  mamat = fill(Inf, flen, size(mus)[1]); # in id order
  valids = vec(sum(mus, dims = 2) .> 0);
  
  if typeof(distances) <: Matrix
  
    fulldist = fill(Inf, size(mus)[1]);
    _fixed_rankmatches!(ranks, fulldist, distances, mus, flen, valids);

  else
    mamat[:, valids] = @views distances[1];
    fbools = sum(mamat .> 0, dims = 2) .> 0;

    _rankmatches!(
      ranks, fbools,
      eachrow(mamat),
    );
  end
  
  return tob
end

function _rankmatches!(ranks, fbools, mamatc)

  # (φ, (fbool, ec, mc)) = collect(enumerate(zip(fbools, mamatc, musc)))[1]

  for (φ, (fbool, ec)) in enumerate(zip(fbools, mamatc))
    if fbool
      # in-place ranking
      # ranks[φ] = StatsBase.ordinalrank(ec);
      # ranks[φ][isinf.(ec)] .= 0; # make zero if not valid

      # really, we want to be able to shut off the mus easily
      # and make it such that we get the mu elements in ranked order
      # by simple indexing:
      # ranks[φ] = @views mc[sortperm(ec)];
  
      # give indices based on ids
      # only give valid ones
      ranks[φ] = sortperm(ec)[1:findfirst(isinf.(ec[sortperm(ec)]))-1]

      # remove impossible mus (for that f)
      # ranks[φ][findfirst(isinf.(ec[sortperm(ec)])):end] .= 0
        # Inf values should never be used, in any model
      # example refinement application:
      # ormus[1+5:end] .= 0
    end
  end
  return ranks
end

function _fixed_rankmatches!(ranks, fulldist, distances, mus, flen, valids)
  for φ in eachindex(1:flen)
    valf = @views mus[:, φ];
    fulldist .= Inf;
    fulldist[valids] = @views distances[:, 1]; # mahalanobis is col 1
    sp = sortperm(fulldist[valf]);
    ranks[φ] = (1:size(mus)[1])[valf][sp];
  end
  return ranks
end
