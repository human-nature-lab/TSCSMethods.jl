# _makegroup_indices.jl

# types for group indices object
STypeMat = SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Base.Slice{Base.OneTo{Int64}}}, false};

STypeMatMis = SubArray{Union{Missing, Float64}, 2, Matrix{Union{Missing, Float64}}, Tuple{Vector{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}

STypeVec = SubArray{Int64, 1, Vector{Int64}, Tuple{Vector{Int64}}, false}
STypeVecBool = SubArray{Bool, 1, Vector{Bool}, Tuple{Vector{Bool}}, false}
# last entry probably depends on original datatype before conversion

"""
    make_groupindices(tvec, treatvec, idvec, uid, fmax, Lmin, cdat)

Get partially overlapping values for treated observations.
"""
function make_groupindices(
  tvec, treatvec, idvec, uid, fmax, Lmin, cdat;
  exvec = nothing
)

  treatvecbool = convert(Vector{Bool}, treatvec);
  tts = sort(unique(tvec[treatvecbool]));
  
  XT = if Missing <: eltype(cdat)
    STypeMatMis
  else
    STypeMat
  end

  tidx = Dict{Tuple{Int, Int}, XT}()

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
      tts, uid, fmax, Lmin,
      tvec, idvec,
      treatvecbool, cdat
    )

    return tidx, ridx, tridx
  else
    _makegroupindices(
      tidx, ridx, tridx,
      tts, uid, fmax, Lmin,
      tvec, idvec,
      treatvecbool, cdat,
      exidx, exvec
    )
    return tidx, ridx, tridx, exidx
  end
  
end


function _makegroupindices(
  tidx, ridx, tridx, tts, uid, fmax, Lmin, tvec, idvec, treatvecbool,
  cdat
)

  # (tt, unit) = collect(Iterators.product(tts, uid))[1]
  @floop for (tt, unit) in Iterators.product(tts, uid)
    @init yesrows = Vector{Bool}(undef, length(tvec))
    getyes!(yesrows, tvec, idvec, tt, fmax, Lmin, unit)

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
  tidx, ridx, tridx, tts, uid, fmax, Lmin, tvec, idvec, treatvecbool,
  cdat, exidx, exvec
)
  @floop for (tt, unit) in Iterators.product(tts, uid)
    @init yesrows = Vector{Bool}(undef, length(tvec))
    getyes!(yesrows, tvec, idvec, tt, fmax, Lmin, unit)

    tidx[(tt, unit)] = @views cdat[yesrows, :];
    ridx[(tt, unit)] = @views tvec[yesrows];
    tridx[(tt, unit)] = @views treatvecbool[yesrows];
    exidx[(tt, unit)] = @views exvec[yesrows];
  end
  return tidx, ridx, tridx, exidx
end

function getyes!(yesrows, tvec, idvec, tt, fmax, Lmin, unit)
  # conditions on which rows of the data to grab
  # iterates over time vector, and picks in-range elements
  # from tt + Lmin to (tt + fmax - 1)
  # lower bound given by covariate matching window minimum,
  # upper bound is given by the upper F requirement on the
  # treatment history matching window
  for k in eachindex(tvec)
    yesrows[k] = ((tvec[k] < tt + fmax) & (tvec[k] >= tt + Lmin)) & (idvec[k] == unit)
  end
  return yesrows
end
