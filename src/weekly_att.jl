"""
helper to sum the bootstrap distributions over some f range
"""
function estsum(bs)
  Z = sum(bs, dims = 2);
  Z = reshape(Z, size(Z)[1])
  return mean(Z), quantile(Z, [0.025, 0.5, 0.975])
end

"""
gives separate sum for each week, respecting days of week, for f window
note that, depending on f window, weeks may have different numbers of days
"""
function weeklyattci(bs, atts, startday::Int64, fmin, fmax)
  # create vector denoting week assignments over cols of bs
  dow = fld.(collect(startday .+ (fmin:fmax)), 7) # from week of primary (0) on
  udow = unique(dow)

  bss = Matrix{Float64}(undef, size(bs)[1], length(udow));
  attsums = Vector{Float64}(undef, length(udow));

  for (i, e) in enumerate(udow)
    lx = findall(dow .== e)
    bss[:, i] = sum(@views(bs[:, lx]), dims = 2)
    attsums[i] = sum(@views(atts[lx]))
  end

  return bss, attsums, udow
end

function processattsum(bss, attsums, udow)
  qtm = Matrix{Float64}(undef, length(attsums), 3)
  for i = eachindex(attsums)
    qtm[i, :] = quantile(bss[:, i], [0.025, 0.5, 0.975])
  end

  sumres = DataFrame(
    week = udow,
    attsum = attsums,
    lwer = qtm[:, 1],
    med = qtm[:, 2],
    uper = qtm[:, 3]
  )
  return sumres
end

function weeklyatt(bs::Matrix{Float64}, res, startday)
  
  res = sort(res, [:f])

  bss, attsums, udow = weeklyattci(
    bs,
    res.att,
    startday,
    minimum(res.f),
    maximum(res.f)
  )
  sumres = processattsum(bss, attsums, udow)
  return sumres
end

function weeklyatt(
  bs::Dict{Int64, Matrix{Float64}},
  res,
  startday # day of the week as Int
)
  res = sort(res, [:stratum, :f])

  SumRes = mkDataFrame(
    Dict(
      :stratum => Int64[],
      :week => Int64[],
      :attsum => Float64[],
      :lwer => Float64[],
      :med => Float64[],
      :uper => Float64[]
    )
  );

  for (k, v) in bs
    resk = @views(res[res.stratum .== k, :])
    bss, attsums, udow = weeklyattci(
      v,
      resk[!, :att],
      startday,
      minimum(resk.f),
      maximum(resk.f)
    )
    sumres = processattsum(bss, attsums, udow)
    sumres[!, :stratum] .= k
    append!(SumRes, sumres)
  end
  return SumRes
end

function plot_att_sum(sumres; savename = "")
  wmin = minimum(sumres.week)

  att_plt = plot(
    sort(sumres, [:week]),
      xintercept = [0],
      y = :week,
      x = :attsum,
      xmin = :lwer, xmax = :uper,
      Geom.point,
      Geom.vline(style = :dot, color = "black", size = [0.2mm]),
      Geom.errorbar,
      Guide.title("ATT (weekly sum)"),
      Guide.ylabel("weeks after the primary"),
      Guide.xlabel("estimate"),
      Coord.Cartesian(ymin = wmin)
  )

  if length(savename) > 0
    draw(PNG(savename, 9inch, 5inch), att_plt)
  end
  return att_plt
end


function plot_att_sum(sumres, stratvar::Symbol; savename = "")
  wmin = minimum(sumres.week)

  att_plt = plot(
      sort(sumres, [:stratum, :week]),
      xintercept = [0],
      y = :week,
      x = :attsum,
      xmin = :lwer, xmax = :uper,
      xgroup = :stratum,
      Guide.title("ATT (weekly sum)" * " by " * string(stratvar)),
      Guide.ylabel("weeks after the primary"),
      Guide.xlabel("estimate"),
      Geom.subplot_grid(
        Geom.point,
        Geom.vline(style = :dot, color = "black", size = [0.2mm]),
        Geom.errorbar,
        Coord.Cartesian(ymin = wmin)
      )
  )

  if length(savename) > 0
    draw(PNG(savename, 9inch, 5inch), att_plt)
  end
  return att_plt
end

function plot_att_sum(sumres, stratvar::Symbol, stratdict; savename = "")
  wmin = minimum(sumres.week)

  sumres.svlabel = getindex.(Ref(stratdict), sumres[!, :stratum])  

  att_plt = plot(
      sort(sumres, [:stratum, :week]),
      xintercept = [0],
      y = :week,
      x = :attsum,
      xmin = :lwer, xmax = :uper,
      xgroup = :svlabel,
      Guide.title("ATT (weekly sum)" * " by " * string(stratvar)),
      Guide.ylabel("weeks after the primary"),
      Guide.xlabel("estimate"),
      Geom.subplot_grid(
        Geom.point,
        Geom.vline(style = :dot, color = "black", size = [0.2mm]),
        Geom.errorbar,
        Coord.Cartesian(ymin = wmin)
      )
  )

  if length(savename) > 0
    draw(PNG(savename, 9inch, 5inch), att_plt)
  end
  return att_plt
end