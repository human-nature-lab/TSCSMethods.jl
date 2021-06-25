# plotting_functions.jl

# matching plots
  
function plot_balance(
  meanbalances, when::String;
  xinch = 89mm,
  yinch = 89mm * 0.75,
  savename = ""
)

  ncv = Symbol("Matching Covariate");

  plt = plot(
    rename(
      sort(meanbalances, [:covariate]),
      :covariate => "Matching Covariate",
    ),
    x = :matchtime,
    y = :meanscore,
    color = ncv,
    Geom.point, Geom.line,
    Guide.xticks(
      ticks = [
        minimum(meanbalances.matchtime):1:maximum(meanbalances.matchtime);
      ],
      orientation = :vertical
    ),
    Guide.title("Covariate Balance (" * when * "-Refinement)"),
    Guide.xlabel("Day Before Treatment"),
    Guide.ylabel("Balance Score"),
    Theme(
      key_position = :right,
      key_label_font_size = 07pt,
      key_title_font_size = 09pt
    )
  )
  if length(savename) > 0
    draw(PNG(savename, xinch, yinch), plt)
  end
  return plt
end

function plot_balance(
  meanbalances,
  stratvar::Symbol,
  when::String;
  savename = "",
  xinch = 89mm * 2,
  yinch = 89mm * 1.5 * 2,
)

  sv = String(stratvar) * " Stratum";

  ncv = Symbol("Matching Covariate");

  plt = plot(
    rename(
      sort(meanbalances, [:covariate, sv]),
      :covariate => "Matching Covariate",
    ),
    x = :matchtime,
    y = :meanscore,
    color = ncv,
    ygroup = sv,
    Guide.title("Covariate Balance " * when * "-Refinement"),
    Guide.xlabel("Day Before Treatment", orientation = :horizontal),
    Guide.ylabel("Balance Score", orientation = :vertical),
    Geom.subplot_grid(
      Geom.point, Geom.line,
      free_y_axis = true,
      Guide.xticks(ticks = [minimum(meanbalances.matchtime):1:maximum(meanbalances.matchtime);])
    ),
    Theme(
      key_position = :right,
      key_label_font_size = 07pt,
      key_title_font_size = 09pt
    )
  )
  if length(savename) > 0
    draw(PNG(savename, xinch, yinch), plt)
  end
  return plt
end

"""
stratdict::Dict points from stratvar to desired label
"""
function plot_balance(
  meanbalances,
  stratvar::Symbol,
  when::String,
  stratdict::Dict;
  savename = "",
  xinch = 89mm * 2,
  yinch = 89mm * 1.5 * 2,
)

  sv = String(stratvar) * " Stratum";

  meanbalances.svlabel = getindex.(Ref(stratdict), meanbalances[!, sv])

  ncv = Symbol("Matching Covariate");

  plt = plot(
    rename(
      sort(meanbalances, [:covariate, sv]),
      :covariate => "Matching Covariate",
    ),
    x = :matchtime,
    y = :meanscore,
    color = ncv,
    ygroup = :svlabel,
    Guide.title("Covariate Balance " * when * "-Refinement"),
    Guide.xlabel("Day Before Treatment"),
    Guide.ylabel("Balance Score"),
    Geom.subplot_grid(
      Geom.point, Geom.line,
      free_y_axis = true,
      Guide.xticks(
        ticks = [
          minimum(meanbalances.matchtime):1:maximum(meanbalances.matchtime);
        ],
        orientation = :vertical
      )
    ),
    Theme(
      key_position = :right,
      key_label_font_size = 07pt,
      key_title_font_size = 09pt
    )
  )
  if length(savename) > 0
    draw(PNG(savename, xinch, yinch), plt)
  end
  return plt
end

function plot_balance(
  meanbalances,
  stratvar::Symbol,
  when::String,
  ygv::Symbol;
  savename = "",
  xinch = 89mm * 2,
  yinch = 89mm * 1.5 * 2,
)

  sv = String(stratvar) * " Stratum";

  ncv = Symbol("Matching Covariate");

  plt = plot(
    rename(
      sort(meanbalances, [:covariate, sv]),
      :covariate => ncv,
    ),
    x = :matchtime,
    y = :meanscore,
    color = ncv,
    ygroup = ygv,
    Guide.title("Covariate Balance " * when * "-Refinement"),
    Guide.xlabel("Day Before Treatment"),
    Guide.ylabel("Balance Score"),
    Geom.subplot_grid(
      Geom.point, Geom.line,
      free_y_axis = true),
      Guide.xticks(ticks = [minimum(meanbalances.matchtime):1:maximum(meanbalances.matchtime);]
      )
    )
  if length(savename) > 0
    draw(PNG(savename, xinch, yinch), plt)
  end
  return plt
end


function plot_balance(
  meanbal_pre,
  meanbal_post;
  savename = ""
)

  meanbal_pre.pre = repeat(["pre-refinement"], nrow(meanbal_pre))
  meanbal_post.pre = repeat(["post-refinement"], nrow(meanbal_post))
  
  # long format for Gadfly
  meanbal = vcat(meanbal_pre, meanbal_post)

  ncv = Symbol("Matching Covariate");

  plt = plot(
    rename(
      sort(meanbalances, [:covariate]),
      :covariate => ncv,
    ),
    meanbal,
    x = :matchtime,
    y = :meanscore,
    color = ncv,
    Scale.y_continuous(minvalue = -0.3, maxvalue = 0.3),
    Guide.title("Covariate Balance"),
    Guide.xlabel("Day Before Treatment"),
    Guide.ylabel("Balance Score"),
    xgroup = :pre,
    Geom.subplot_grid(
      Geom.point,
      Geom.line,
      Guide.xticks(ticks = [minimum(meanbalances.matchtime):1:maximum(meanbalances.matchtime);])
    )
  )

  if length(savename) > 0
    draw(PNG(savename, 18inch, 5inch), plt)
  end
  return plt
end

# estimation plots

"""
    plot_att(atts, estimator; savename = "")

"""
function plot_att(
  atts, estimator;
  xinch = 89mm,
  yinch = 89mm * 0.75,
  savename = "")
  fmin = minimum(atts.f)
  fmax = maximum(atts.f)

  att_plt = plot(
    atts,
    layer(
      x = :f,
      y = :att,
      ymin = :lwer,
      ymax = :uper,
      Geom.point,
      Geom.errorbar
    ),
    layer(
      yintercept = [0],
      Geom.hline(style = :dot, color = "black", size = [0.5mm])
    ),
    Guide.xlabel("Day After Treatment"),
    Guide.ylabel("Estimate"),
    Guide.title(estimator),
    Guide.xticks(ticks = collect(fmin : fmax), orientation = :vertical)
  )

  if length(savename) > 0
    draw(PNG(savename, xinch, yinch), att_plt)
  end
  return att_plt
end

"""
    plot_att(
      atts,
      stratvar::Symbol,
      savename = "",
      xinch = 15inch,
      yinch = 6inch,
      treatment = :treatment
    )

Plot restricted estimates, without key for strata.

"""
function plot_att(
  atts,
  stratvar::Symbol,
  estimator;
  savename = "",
  xinch = 89mm * 2,
  yinch = 89mm * 0.75 * 1.5 * 2,
  treatment = :treatment
)

  # sv = String(stratvar) * "_stratum";

  ttl = estimator * " by " * String(replace(String(stratvar), "_" => " "))

  fmin = minimum(atts.f)
  fmax = maximum(atts.f)

  plt = plot(
    sort(atts, [:stratum, :f]),
    xintercept = [0],
    y = :f,
    x = :att,
    xmin = :lwer, xmax = :uper,
    ygroup = :stratum,
    Guide.title(ttl),
    Guide.xlabel("Day After Treatment"),
    Guide.ylabel("Estimate"),
    Geom.subplot_grid(
      Geom.vline(style = :dot, color = "black", size = [0.5mm]),
      Geom.point,
      Geom.errorbar,
      free_y_axis = true,
      # free_x_axis = true,
      Guide.yticks(ticks = fmin : 1 : fmax)
    )
  )
  if length(savename) > 0
    draw(PNG(savename, xinch, yinch), plt)
  end
  return plt
end

"""
    plot_att(
      atts,
      stratvar::Symbol,
      stratdict::Dict;
      savename = "",
      xinch = 15inch,
      yinch = 6inch,
      treatment = :treatment
    )

Plot restricted estimates, with key for strata.

"""
function plot_att(
  atts,
  stratvar::Symbol,
  stratdict::Dict,
  estimator;
  savename = "",
  xinch = 89mm * 2,
  yinch = 89mm * 0.75 * 1.5 * 2,
  treatment = :treatment
)

  ttl = estimator * " by " * String(replace(String(stratvar), "_" => " "))

  fmin = minimum(atts.f)
  fmax = maximum(atts.f)

  # sv = String(stratvar) * "_stratum";
  atts.svlabel = getindex.(Ref(stratdict), atts[!, :stratum])

  plt = plot(
    sort(atts, [:stratum, :f]),
    yintercept = [0],
    x = :f,
    y = :att,
    ymin = :lwer,
    ymax = :uper,
    ygroup = :svlabel,
    Guide.title(ttl),
    Guide.xlabel("Day After Treatment"),
    Guide.ylabel("Estimate"),
    Geom.subplot_grid(
      Geom.point,
      Geom.hline(style = :dot, color = "black", size = [0.5mm]),
      Geom.errorbar,
      free_y_axis = true,
      Guide.xticks(
        ticks = fmin : 1 : fmax,
        orientation = :vertical
      )
    )
  )
  if length(savename) > 0
    draw(PNG(savename, xinch, yinch), plt)
  end
  return plt
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