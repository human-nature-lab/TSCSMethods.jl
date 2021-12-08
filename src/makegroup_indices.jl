# _makegroup_indices.jl

function _makegroupindices(
  tidx, ridx, tridx, tts, uid, fmin, fmax, mmin, tvec, idvec, treatvec, cdat
)
  @floop for (tt, unit) in Iterators.product(tts, uid)
    @init yesrows = Vector{Bool}(undef, length(tvec))
    getyes!(yesrows, tvec, idvec, tt, fmin, fmax, mmin, unit)

    tidx[(tt, unit)] = @views cdat[yesrows, :];
    ridx[(tt, unit)] = @views tvec[yesrows];
    tridx[(tt, unit)] = @views treatvec[yesrows];
  end
  return tidx, ridx, tridx
end

function getyes!(yesrows, tvec, idvec, tt, fmin, fmax, mmin, unit)
  for k in eachindex(tvec)
    yesrows[k] = ((tvec[k] < tt + fmax) & (tvec[k] >= tt + fmin + mmin)) & (idvec[k] == unit)
  end
  return yesrows
end

"""
    make_groupindices(tvec, treatvec, idvec, uid, fmin, fmax, mmin, cdat)

Get partially overlapping values for treated observations.
"""
function make_groupindices(tvec, treatvec, idvec, uid, fmin, fmax, mmin, cdat)
  #  14.354192 seconds
  # (3.11 M allocations: 283.229 MiB, 14.99% gc time, 7.10% compilation time)
  STypeMat = SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Base.Slice{Base.OneTo{Int64}}}, false};
  STypeVec = SubArray{Int64, 1, Vector{Int64}, Tuple{Vector{Int64}}, false}

  tts = sort(unique(tvec[treatvec .== 1]));
  
  tidx = Dict{Tuple{Int, Int}, STypeMat}();
  ridx = Dict{Tuple{Int, Int}, STypeVec}();
  tridx = Dict{Tuple{Int, Int}, STypeVec}();

  sizehint!(tidx, length(Iterators.product(tts, uid)));
  sizehint!(ridx, length(Iterators.product(tts, uid)));
  sizehint!(tridx, length(Iterators.product(tts, uid)));

  _makegroupindices(
    tidx, ridx, tridx, tts, uid, fmin, fmax, mmin, tvec, idvec, treatvec, cdat
  )
  
  return tidx, ridx, tridx
end
