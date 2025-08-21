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
  @unpack mus, distances, ranks = tob;

  # Pre-compute validity mask once
  valids = vec(sum(mus, dims = 2) .> 0);
  
  if typeof(distances) <: Matrix
    # Reuse single buffer instead of creating new one each time
    fulldist = fill(Inf, size(mus)[1]);
    _fixed_rankmatches!(ranks, fulldist, distances, mus, flen, valids);
  else
    # Streaming approach: avoid full matrix allocation
    _streaming_rankmatches!(ranks, distances, mus, flen, valids);
  end
  
  return tob
end

function _rankmatches!(ranks, fbools, mamatc)

  # (φ, (fbool, ec, mc)) = collect(enumerate(zip(fbools, mamatc, musc)))[1]

  for (φ, (fbool, ec)) in enumerate(zip(fbools, mamatc))
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
        ranks[φ] = sp[1:last_finite_idx]
      end
    end
  end
  return ranks
end

function _streaming_rankmatches!(ranks, distances, mus, flen, valids)
  # Streaming approach: process one f at a time without full matrix
  for φ in 1:flen
    valf = @views mus[:, φ]
    if any(valf)
      # Get distances for this f, only for valid matches
      valid_and_allowed = valids .& valf
      if any(valid_and_allowed)
        dist_f = @views distances[1][:, 1][valid_and_allowed]  # mahalanobis col 1
        valid_indices = findall(valid_and_allowed)
        
        # Single-pass ranking: sort indices by distance
        sorted_indices = valid_indices[sortperm(dist_f)]
        ranks[φ] = sorted_indices
      end
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
