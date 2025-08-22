# ranking.jl

## ranking for the sliding window

function rank!(matches, flen)
  Threads.@threads :greedy for i in 1:length(matches)
    tob = @views matches[i];
    rankmatches!(tob, flen)
  end
  return matches
end

function rankmatches!(tob, flen)
  (; eligible_matches, distances, match_rankings) = tob;

  # Pre-compute validity mask once
  valids = vec(sum(eligible_matches, dims = 2) .> 0);
  
  if typeof(distances) <: Matrix
    # Reuse single buffer instead of creating new one each time
    fulldist = fill(Inf, size(eligible_matches)[1]);
    _fixed_rankmatches!(match_rankings, fulldist, distances, eligible_matches, flen, valids);
  else
    # Streaming approach: avoid full matrix allocation
    _streaming_rankmatches!(match_rankings, distances, eligible_matches, flen, valids);
  end
  
  return tob
end

function _rankmatches!(ranks, fbools, mamatc)

  # (window_index, (fbool, ec, mc)) = collect(enumerate(zip(fbools, mamatc, musc)))[1]

  for (window_index, (fbool, ec)) in enumerate(zip(fbools, mamatc))
    if fbool
      # Single-pass ranking: sort indices by distance, stop at first Inf
      sp = sortperm(ec)
      
      # Find first infinite value efficiently
      last_finite_idx = 0
      for i in 1:length(sp)
        if isinf(ec[sp[i]])
          break
        end
        last_finite_idx = i
      end
      
      # Only keep finite distances
      if last_finite_idx > 0
        ranks[window_index] = sp[1:last_finite_idx]
      end
    end
  end
  return ranks
end

function _streaming_rankmatches!(match_rankings, distances, eligible_matches, flen, valids)
  # Streaming approach: process one f at a time without full matrix
  for window_index in 1:flen
    valf = @views eligible_matches[:, window_index]
    if any(valf)
      # Get distances for this f, only for valid matches
      valid_and_allowed = valids .& valf
      if any(valid_and_allowed)
        dist_f = @views distances[1][:, 1][valid_and_allowed]  # mahalanobis col 1
        valid_indices = findall(valid_and_allowed)
        
        # Single-pass ranking: sort indices by distance
        sorted_indices = valid_indices[sortperm(dist_f)]
        match_rankings[window_index] = sorted_indices
      end
    end
  end
  return match_rankings
end

function _fixed_rankmatches!(match_rankings, fulldist, distances, eligible_matches, flen, valids)
  for window_index in eachindex(1:flen)
    valf = @views eligible_matches[:, window_index];
    fulldist .= Inf;
    fulldist[valids] = @views distances[:, 1]; # mahalanobis is col 1
    sp = sortperm(fulldist[valf]);
    match_rankings[window_index] = (1:size(eligible_matches)[1])[valf][sp];
  end
  return match_rankings
end
