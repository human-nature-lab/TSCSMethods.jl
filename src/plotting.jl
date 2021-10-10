# from Gadfly documentation
# <http://gadflyjl.org/stable/gallery/scales/#[Scale.color_discrete_hue](@ref)>
function gen_colors(n)
    cs = Colors.distinguishable_colors(
      n,
      [colorant"#FE4365", colorant"#eca25c"],
      lchoices = Float64[58, 45, 72.5, 90],
      transform = c -> deuteranopic(c, 0.1),
      cchoices = Float64[20,40],
      hchoices = [75,51,35,120,180,210,270,310]
    )
    return convert(Vector{Color}, cs)
  end

function mk_covpal(variables::Vector{Symbol})
  varcol = Dict{Symbol, RGB}();

  pal = gen_colors(length(variables));

  for (i, variable) in enumerate(variables)
    varcol[variable] = pal[i]
  end
  return varcol
end

function mkmname(ttl, cspec, outc, cal; sv = "")

  if isnothing(sv)
    sv = ""
  end

  calstr = cal == true ? "cal" : "";
  outcstr = string(outc);
  
  mname = ttl * cspec * "_" * outcstr * "_" * sv * "_" * calstr * "_model";

  return mname
end

"""
ax_att()

Generate and fill in the axis for the ATT plot. Used in conjunction with model_pl().
  
There is no standalone att_plot().
"""
function ax_att(
  fsub,
  atts;
  attl = "",
  outcome = :death_rte,
  CIbnds = [Symbol("2.5%"), Symbol("97.5%")]
)

  # paper specific
  if outcome == :death_rte
    xt = collect(10:3:40)
    ocolor = gen_colors(3)[3]
    olab = "Death Rate"
  elseif outcome == :case_rte
    xt = collect(3:3:40)
    ocolor = gen_colors(5)[5]
    olab = "Case Rate"
  # general
  else
    fmin = minimum(atts.f); fmax = maximum(atts.f)
    xt = collect(fmin:(fmax-fmin)/10:fmax)
    ocolor = gen_colors(3)[3]
    olab = string(outcome)
  end

  ax = Axis(
    fsub,
    title = attl,
    xticks = xt,
    xlabel = "Day",
    ylabel = "ATT Estimate",
    xminorgridvisible = true,
    xminorticksvisible = true,
    xminorticks = IntervalsBetween(3)
  )

  rb = rangebars!(
    atts.f,
    atts[!, CIbnds[1]], atts[!, CIbnds[2]],
    color = ocolor;
    whiskerwidth = 10,
    label = olab
  )

  # plot position scatters so low and high errors can be discriminated
  scatter!(atts.f, atts.att, markersize = 5, color = :black)

  return ax, rb
end

"""
    makeseries(cbi; variablecolors::Union{Dict, Nothing} = nothing)

Convert a grandbalances dictionary into a series, with labels and colors, to plot. Output goes into series!(), inputs a grandbalances object for a non-stratified analysis, or a single stratum.

variablecolors allows input of a custom color set.
"""
function makeseries(cbi; variablecolors::Union{Dict{Symbol, RGB}, Nothing} = nothing)
  # length of time-varying coeff
  clen = maximum([length(v) for v in values(cbi)]);
  rlen = length(keys(cbi));
  labs = sort([k for k in keys(cbi)]);
  
  servals = Matrix{Float64}(undef, rlen, clen);
  for (r, lab) in enumerate(labs)
    vals = cbi[lab]
    if length(vals) == 1
      servals[r, :] .= vals;
    else
      servals[r, :] = vals
    end
  end

  if !isnothing(variablecolors)
    varcol = mk_covpal(variablecolors)
  else
    varcol = variablecolors
  end
  
  sercols = [varcol[lab] for lab in labs];
  
  return servals, string.(labs), sercols
end

function ax_cb(
  fsub,
  cbi,
  variablecolors;
  cbttl = "",
  step = 5
)

  servals, serlabs, sercols = makeseries(cbi; variablecolors = variablecolors);

  axcb = Axis(
    fsub,
    title = cbttl,
    xticks = collect(range(-50, stop = 0; step = step)), # generalize
    xlabel = "Day",
    ylabel = "Balance Score",
    xminorgridvisible = true,
    xminorticksvisible = true,
    xminorticks = IntervalsBetween(step)
  );

  ser = series!(
    axcb,
    collect(range(-50, -1; step =1)), # generalize
    servals,
    labels = serlabs,
    markersize = 5,
    color = sercols
  )

  return axcb, ser
end

function load_mod(
  ttl,
  epth,
  outcome,
  cspec,
  cal;
  stratvar = nothing
)

  jext = ".jld2";

  mname = mkmname(ttl, cspec, outcome, cal; sv = stratvar);
  model = JLD2.load_object(epth * mname * jext);

  println(mname)

  return model
end

function plot_cbs(
  model1::Union{cicmodel, calipercicmodel},
  model2::refinedcicmodel;
  labels = Dict{Int, String}(),
  variablecolors = nothing,
  fw = 700,
  fl = 300,
  spath = nothing
)

  txt1 = "Pre-Refinement Balance";
  txt2 = "Post-Refinement Balance";

  f1 = []
  f2 = []

  f1 = plot_cb(
    model1;
    labels = labels,
    variablecolors = variablecolors,
    fw = fw, fl = fl,
    spath = spath
  )

  supertitle = f1[0, :] = Label(
    f1,
    txt1,
    color = (:black, 0.25)
  )

  f2 = plot_cb(
    model2;
    labels = labels,
    variablecolors = variablecolors,
    fw = fw, fl = fl,
    spath = spath
  )

  supertitle = f2[0, :] = Label(
    f2,
    txt2,
    color = (:black, 0.25)
  )
  
  return [f1, f2]
end

"""
    plot_cb()

Construct a pre- or post-refinement balance plot. Handles two cases:
1) single plot, not stratified
2) one plot for each stratum
"""
function plot_cb(
  model;
  labels = Dict{Int, String}(),
  variablecolors = nothing,
  fw = 700, fl = 300,
  spath = nothing
)

  strat = model.stratifier != Symbol("")
  svlab = string(model.stratifier)

  if !strat
    lf = 1; wf = 1
    ns = 1
  elseif strat
    pdict = make_pdict();

    # sl = svlab * " Stratum";
    sn = sort(unique(
      sort([k for k in keys(model.grandbalances)])
    ));

    lf, wf, ns = fdims(sn)
  end

  f = Figure(resolution = (wf * fw, lf * fl + fl));

  if strat

    if isempty(labels)
      error("supply strata labels")
    end
    
    for (i, s) in enumerate(sn)

      cbi = model.grandbalances[s];

      fpos = pdict[i];
      axc, ser = ax_cb(f[fpos...][1, 1], cbi, variablecolors);

      label_tt = Label(
        f,
        labels[s],
        halign = :center)
      f[fpos..., Top()] = label_tt
    end
  else

    fpos = [1, 1]
    axc, ser = ax_cb(f[fpos...][1,1], cb, variablecolors);

    hm_sublayout = GridLayout()
    f[1, 1] = hm_sublayout

  end

  nbanks = ns == 1 ? 2 : 4;

  lcb = Legend(f, axc, "Covariates", framevisible = false, nbanks = nbanks)

  hm_sublayout2 = GridLayout()
  f[lf + 1, :] = hm_sublayout2

  hm_sublayout2[:v] = [lcb]
  # colsize!(hm_sublayout2, 1, Aspect(1,1))

  if !isnothing(spath)

    cal = !isempty(model.caliper)
    mname = mkmname(model.title, "", model.outcome, cal; sv = svlab)
    when = typeof(refinedcicmodel) ? "Post" : "Pre";
    
    save(
      spath * "model_" * when * "_" * mname * ".png",
      f
    )
  end

  return f
end

pl_ratio(mo) = mo == :death_rate ? Relative(13/31) : Relative(13/38)

function model_pl(
  model::AbstractCICModel;
  labels = Dict{Int, String}(),
  variablecolors = nothing,
  fw = 700, fl = 300
)

  crel = pl_ratio(model.outcome);

  # load model function
  # from
  # ttl = "voting";
  # epth = "../covid-19-voting/outp_voting_deaths/";

  strat = !(model.stratifier == Symbol(""))
  # svlabel = string(model.stratifier)
  # svs = svlabel * " Stratum";


  if !strat
    lf = 1; wf = 1
    ns = 1
  elseif strat
    pdict = make_pdict()

    # sl = svlabel * " Stratum";
    sn = sort(unique(
      sort([k for k in keys(model.grandbalances)])
    ));

    lf, wf, ns = fdims(sn)
  end

  f = Figure(resolution = (wf * fw, lf * fl + fl));

  if strat

    if isempty(labels)
      error("supply strata labels")
    end
    
    for (i, s) in enumerate(sn)

      cbi = model.grandbalances[s];
      resi = model.results[model.results.stratum .== s, :];
      
      fpos = pdict[i];
      axa, rb = ax_att(f[fpos...][1,1], resi; outcome = model.outcome);
      axc, ser = ax_cb(f[fpos...][1,2], cbi, variablecolors; step = 10);

      axs = [axc, axa];

      hm_sublayout = GridLayout()
      f[fpos...][1, 1:2] = hm_sublayout

      hm_sublayout[:h] = axs
      colsize!(hm_sublayout, 1, crel)

      label_tt = Label(
        f,
        labels[s],
        halign = :center)
      f[fpos..., Top()] = label_tt
    end
  else
    fpos = [1,1]
    axa, rb = ax_att(f[fpos...][1,1], model.results; outcome = model.outcome);
    axc, ser = ax_cb(f[fpos...][1,2], model.grandbalances, variablecolors);

    axs = [axc, axa];

    hm_sublayout = GridLayout()
    f[1, 1:2] = hm_sublayout

    hm_sublayout[:h] = axs
    colsize!(hm_sublayout, 1, crel)
  end

  nbanks = ns == 1 ? 2 : 4;

  lcb = Legend(f, axc, "Covariates", framevisible = false, nbanks = nbanks)
  latt = Legend(f, axa, "Outcome", framevisible = false)

  hm_sublayout2 = GridLayout()
  f[lf + 1, :] = hm_sublayout2

  hm_sublayout2[:v] = [lcb, latt]
  # colsize!(hm_sublayout2, 1, Aspect(1,1))

  return f
end

function fdims(sn)
  ns = length(sn);
  if ns == 1
    lf = 1; wf = 1
  elseif ns == 2
    lf = 1; wf = 2
  elseif (ns == 3) | (ns == 4)
    lf = 2; wf = 2
  elseif ns == 5
    lf = 3; wf = 2
  elseif ns == 8
    lf = 4; wf = 2
  end
  return lf, wf, ns
end

"""
    make_pdict()

Manually-specified figure placement, for up to eight strata.
"""
function make_pdict()

  pdict = Dict{Int64, Vector{Int64}}();
  pdict[1] = [1,1]
  pdict[2] = [1,2]
  pdict[3] = [2,1]
  pdict[4] = [2,2]
  pdict[5] = [3,1]
  pdict[6] = [3,2]
  pdict[7] = [4,1]
  pdict[8] = [4,2]

  return pdict
end

"""
    plot_modelset(model_path; variablecolors = nothing, base_savepath = "")

Generate the plots, in a new directory, for a set of models in some model set file. Base_savepath should end in "/".
"""
function plot_modelset(
  model_path;
  variablecolors = nothing,
  base_savepath = "" # ends in /
)

  mpnames = [
    "model_plot.png"
    "refined_plot.png"
    "caliper_plot.png"
    "refinedcaliper_plot.png"
  ];

  # cc, ccr, cal, calr, labels
  objet = load_object(model_path); # better to load outside?
  labels = objet[length[objet]]

  dirn = base_savepath * name_model(cc);
  mkdir(dirn);

  for i in 1:length(objet)-1
    mod = objet[i]
    if !isnothing(mod)
        mp = model_pl(
        mod;
        labels = labels,
        variablecolors = variablecolors
      )
      
      save(
        base_savepath * dirn * "/" * mpnames[i],
        mp
      )
    end
  end
end