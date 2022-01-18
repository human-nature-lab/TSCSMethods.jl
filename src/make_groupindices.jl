# _makegroup_indices.jl

function _makegroupindices(
  tidx, ridx, tridx, tts, uid, fmin, fmax, mmin, tvec, idvec, treatvecbool, cdat
)
  @floop for (tt, unit) in Iterators.product(tts, uid)
    @init yesrows = Vector{Bool}(undef, length(tvec))
    getyes!(yesrows, tvec, idvec, tt, fmin, fmax, mmin, unit)

    tidx[(tt, unit)] = @views cdat[yesrows, :];
    ridx[(tt, unit)] = @views tvec[yesrows];
    tridx[(tt, unit)] = @views treatvecbool[yesrows];
  end
  return tidx, ridx, tridx
end

"""
version with exposure
"""
function _makegroupindices(
  tidx, ridx, tridx, tts, uid, fmin, fmax, mmin, tvec, idvec, treatvecbool,
  cdat, exidx, exvec
)
  @floop for (tt, unit) in Iterators.product(tts, uid)
    @init yesrows = Vector{Bool}(undef, length(tvec))
    getyes!(yesrows, tvec, idvec, tt, fmin, fmax, mmin, unit)

    tidx[(tt, unit)] = @views cdat[yesrows, :];
    ridx[(tt, unit)] = @views tvec[yesrows];
    tridx[(tt, unit)] = @views treatvecbool[yesrows];
    exidx[(tt, unit)] = @views exvec[yesrows];
  end
  return tidx, ridx, tridx, exidx
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
function make_groupindices(
  tvec, treatvec, idvec, uid, fmin, fmax, mmin, cdat;
  exvec = nothing
)
  #  14.354192 seconds
  # (3.11 M allocations: 283.229 MiB, 14.99% gc time, 7.10% compilation time)

  STypeMat = SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Base.Slice{Base.OneTo{Int64}}}, false};
  STypeVec = SubArray{Int64, 1, Vector{Int64}, Tuple{Vector{Int64}}, false}
  STypeVecBool = SubArray{Bool, 1, Vector{Bool}, Tuple{Vector{eltype(treatvec)}}, false}
  # last entry probably depends on original datatype before conversion

  treatvecbool = convert(Vector{Bool}, treatvec);
  tts = sort(unique(tvec[treatvecbool]));
  
  tidx = Dict{Tuple{Int, Int}, STypeMat}();
  ridx = Dict{Tuple{Int, Int}, STypeVec}();
  tridx = Dict{Tuple{Int, Int}, STypeVecBool}();
  exidx = Dict{Tuple{Int, Int}, STypeVec}();

  sizehint!(tidx, length(Iterators.product(tts, uid)));
  sizehint!(ridx, length(Iterators.product(tts, uid)));
  sizehint!(tridx, length(Iterators.product(tts, uid)));
  sizehint!(exidx, length(Iterators.product(tts, uid)));

  if isnothing(exvec)
    _makegroupindices(
      tidx, ridx, tridx,
      tts, uid, fmin, fmax, mmin,
      tvec, idvec,
      treatvecbool, cdat
    )

    return tidx, ridx, tridx
  else
    _makegroupindices(
      tidx, ridx, tridx,
      tts, uid, fmin, fmax, mmin,
      tvec, idvec,
      treatvecbool, cdat,
      exidx, exvec
    )
    return tidx, ridx, tridx, exidx
  end
  
end