# plotting_functions.jl

using Gadfly
import Cairo, Fontconfig

# matching plots

function plot_balance(meanbalances)
  return plot(
    meanbalances,
    x = :matchtime,
    y = :meanscore,
    color = :covariate,
    Geom.point, Geom.line,
    #Guide.yticks(ticks = [-0.3:0.1:0.3;]),
    Guide.title("Covariate Balance"),
    Guide.xlabel("Match Period"),
    Guide.ylabel("Standardized Balance Score")
  )
end
  
function plot_balance(meanbalances, savename::String)
  plt = plot(
    meanbalances,
    x = :matchtime,
    y = :meanscore,
    color = :covariate,
    Geom.point, Geom.line,
    #Guide.yticks(ticks = [-0.3:0.1:0.3;]),
    Guide.title("Covariate Balance"),
    Guide.xlabel("Match Period"),
    Guide.ylabel("Standardized Balance Score")
  )
  draw(PNG(savename, 9inch, 5inch), plt)
  return plt
end

function plot_balance(meanbalances, stratvar::Symbol, savename::String)
  sv = String(stratvar) * "_stratum";
  plt = plot(
    meanbalances,
    x = :matchtime,
    y = :meanscore,
    color = :covariate,
    ygroup = sv,
    #Guide.yticks(ticks = [-0.3:0.1:0.3;]),
    Guide.title("Covariate Balance"),
    Guide.xlabel("Match Period"),
    Guide.ylabel("Standardized Balance Score"),
    Geom.subplot_grid(Geom.point, Geom.line)
  )
  draw(PNG(savename, 10inch, 18inch), plt)
  return plt
end
  
function plot_balance(meanbal_pre, meanbal_post; savename = "")
  meanbal_pre.pre = repeat(["pre-refinement"], nrow(meanbal_pre))
  meanbal_post.pre = repeat(["post-refinement"], nrow(meanbal_post))
  
  # long format for Gadfly
  meanbal = vcat(meanbal_pre, meanbal_post)

  plt = plot(
    meanbal,
    x = :matchtime,
    y = :meanscore,
    color = :covariate,
    Scale.y_continuous(minvalue = -0.3, maxvalue = 0.3),
    Guide.title("Covariate Balance"),
    Guide.xlabel("Match Period"),
    Guide.ylabel("Standardized Balance Score"),
    xgroup = :pre,
    Geom.subplot_grid(Geom.point, Geom.line)
  )

  if length(savename) > 0
    draw(PNG(savename, 18inch, 5inch), plt)
  end
  return plt
end

# estimation plots

function plot_att(atts, fmin, savename)
  fmin = minimum(atts.f)

  att_plt = plot(
      atts,
      x = :f, y = :att,
      ymin = :lwer, ymax = :uper,
      Geom.point,
      Geom.errorbar,
      Guide.title("avg. effect of treatment on the treated"),
      Guide.xlabel("f"),
      Guide.ylabel("estimate"),
      Coord.Cartesian(xmin=fmin)
  )

  draw(PNG(savename, 9inch, 5inch), att_plt)
  return att_plt
end

function plot_att(atts, fmin, savename, groupkey, groupvar; dirloc = "")

  fmin = minimum(atts.f)

  atts = leftjoin(atts, groupkey, on = :stratum => :cartgroup)

  strata = unique(atts[groupvar])

  for si in strata

    tt = "avg. effect of treatment on the treated " * string(si)

    attsi = atts[findall(atts[groupvar] .== si), :]

    att_plt = plot(
    attsi,
    x = :f, y = :att,
    ymin = :lwer, ymax = :uper,
    Geom.point,
    Geom.errorbar,
    Guide.title(tt),
    Guide.xlabel("f"),
    Guide.ylabel("estimate"),
    Coord.Cartesian(xmin=fmin)
    )

    sn = dirloc * string(si) * " " * savename

    draw(PNG(sn, 9inch, 5inch), att_plt)
  end
  return att_plt
end

function plot_att_strat(atts, fmin, savename, dirloc)

  S = unique(atts.stratum);
  for s in S

    tt = "avg. effect of treatment on the treated " * string(s)

    attsi = atts[findall(atts.stratum .== s), :]

    att_plt = plot(
      attsi,
      x = :f, y = :att,
      ymin = :lwer, ymax = :uper,
      Geom.point,
      Geom.errorbar,
      Guide.title(tt),
      Guide.xlabel("f"),
      Guide.ylabel("estimate"),
      Coord.Cartesian(xmin=fmin)
      )

    sn = dirloc * string(s) * " " * savename

    draw(PNG(sn, 9inch, 5inch), att_plt)
  end
end

"""
    plot_mdistances(matches5, calvars; savenme = "")

histogram plot the maha distances: overall, and for each identified
caliper variable

if savenme is left empty, plot will display but not save to file
"""
function plot_mdistances(matches5, calvars; savenme = "")
  poss5 = @views(matches5[matches5.possible .== 1, :]);
  calvarsnmes = [Symbol(String(calvar) * "_mdist") for calvar in calvars];
  poss5 = stack(poss5, vcat(:mdist, calvarsnmes));
  mdisthist = plot(
    poss5,
    x = :value,
    xgroup = :variable,
    Geom.subplot_grid(Geom.histogram),
    Guide.title("Mahalanabois Dists: treated obs to possible match"),
    );
  
  if length(savenme) > 0
    draw(PNG(savenme, 9inch, 5inch), mdisthist)
  end
  return mdisthist
end
