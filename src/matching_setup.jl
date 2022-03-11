# matching_setup.jl

function observe(datt, datid, dattrt)
  obslen = sum(dattrt)
  v = Vector{Tuple{Int, Int}}(undef, obslen);
  cnt = 0
  for (τ, unit, z) in zip(datt, datid, dattrt)
    if z > 0
      cnt += 1
      v[cnt] = (τ, unit)
    end
  end

  return sort(v), unique(datid)
end

function make_matches(obslen, idlen, flen)
  matches = Vector{Tob}(undef, obslen);
  _make_matches!(matches, obslen, idlen, flen);
  return matches
end

function _make_matches!(matches, obslen, idlen, flen; sliding = false)
  for i in 1:obslen
    matches[i] = Tob(
      mus = fill(true, idlen, flen), # DEFAULT TRUE -- ACTIVELY REMOVE MATCHES
      fs = Vector{Vector{Bool}}(undef, idlen),
      ranks = Dict{Int, Vector{Int}}()
    );
  end;
end
