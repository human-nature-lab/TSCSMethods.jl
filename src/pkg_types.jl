using Parameters

#=
model object
* the data is kept separately
=#
ReforDataFrame = Union{Base.RefValue{DataFrame}, DataFrame}; # type for model
StratDict = Dict{Int64, Dict{Symbol, Float64}};

@with_kw mutable struct cicmodel
  title::String = "standard"
  id::Symbol
  t::Symbol
  outcome::Symbol
  treatment::Symbol
  matchingcovar::Vector{Symbol}
  caliper::Union{Dict{Symbol, Float64}, StratDict} = Dict{Symbol, Float64}()
  fmin::Int64
  fmax::Int64
  stratvar::Symbol = Symbol()
  tpoint::Int64
  tmin::Int64 # gives length of matching period with tpoint
  refinementnum::Int64 = Int64(0)
  # store either dataframe or reference to dataframe, dep on use cases
  matches::ReforDataFrame = Base.RefValue{DataFrame}()
  matches5::ReforDataFrame = Base.RefValue{DataFrame}()
  data::ReforDataFrame = Base.RefValue{DataFrame}()
  # balances pre refinement
  balances_pre::DataFrame = DataFrame()
  meanbalances_pre::DataFrame = DataFrame()
  # post refinement
  balances_post::DataFrame = DataFrame()
  meanbalances_post::DataFrame = DataFrame()
  # one results per model, so store the DF
  boot_iterations::Int64 = 500
  results::DataFrame = DataFrame() # post refinement results
  results_pre::DataFrame = DataFrame();
  boot_estimates::Union{Array{Float64}, Dict{Int64,Array{Float64,2}}} = Array{Float64}(undef)
    # need to alter estimation output to make this happen
  # plot slots
  # covariate balance pre and post refinement
  pl_cb_pre::Plot = Plot()
  pl_cb_post::Plot = Plot()
  # att estimates pre and post refinement
  pl_att_pre::Plot = Plot()
  pl_att_post::Plot = Plot()
  # total number of treated observations
  treatednum::Union{Int64, Dict{Int64, Int64}} = Int64(0)
  # number of treated left over after filtering or caliper
  treatedleft::Union{Int64, Dict{Int64, Int64}} = Int64(0)
end


# function wrappers for type

function matching!(model::cicmodel, variancesonly::Bool;
  refine = true
)
  
  model.matches = matching(
    model.data[], model.matchingcovar,
    model.id, model.t,
    model.fmin, model.fmax,
    model.treatment,
    model.tpoint,
    model.matchingcovar, # use these for caliper distance calc
    variancesonly);

  model = refine!(model);
  
  return model
end

"""
for separate refinement on a given match set
this avoids time-consuming recalculation of the possible match set
"""
function refine!(model::cicmodel)
  model.matches5 = refine(model.matches, model.refinementnum
);

  actualcal = copy(model.caliper)

  if model.stratvar != Symbol("")
    model.caliper = StratDict()
    caliper!(model) # to fill in treated unit nums, using empty caliper
    model.caliper = copy(actualcal)
  elseif model.stratvar == Symbol("")
    model.caliper = Dict{Symbol, Float64}()
    caliper!(model) # to fill in treated unit nums, using empty caliper
    model.caliper = copy(actualcal)
  end
  
  return model
end

function getbalance!(model::cicmodel; post = true)

  if post == true
    b = :balances_post
    mb = :meanbalances_post
    wmatches = :matches5
  elseif post == false
    b = :balances_pre
    mb = :meanbalances_pre
    wmatches = :matches
  end

  c1 = model.stratvar == Symbol("")
  # c2 = typeof(model.caliper) .== Dict{Symbol,Float64}
  
  if c1

    a1, a2 = getbalance(
      getfield(model, wmatches),
      model.matchingcovar,
      model.data[],
      model.tmin,
      model.tpoint,
      model.id,
      model.t,
      model.treatment
    );
  
  elseif !c1

    a1, a2 = getbalance_restricted(
      getfield(model, wmatches),
      model.stratvar,
      model.matchingcovar,
      model.data[],
      model.tmin,
      model.tpoint,
      model.id,
      model.t,
      model.treatment
    );

  end

  setfield!(model, b, a1)
  setfield!(model, mb, a2)

  return model
end

function estimate!(
  model::cicmodel;
  post::Bool = true,
  ATT::Bool = true,
  outboot::Bool = false
)

  if post == true
    wmatches = :matches5
    when = :results
  else
    wmatches = :matches
    when = :results_pre
  end

  c1 = model.stratvar == Symbol("")

  if c1
    
    res = standard_estimation(
      model.boot_iterations,
      getfield(model, wmatches),
      model.fmin:model.fmax, # nothing fancier will work yet
      model.tpoint,
      model.data[],
      model.id, model.t, model.outcome,
      ATT,
      outboot
    );
  elseif !c1
    res = restricted_estimation(
      model.boot_iterations,
      getfield(model, wmatches),
      model.fmin:model.fmax,
      model.tpoint,
      model.data[],
      model.stratvar,
      id, t, outcome,
      ATT,
      outboot
    );
  end
  
  if outboot == false
    setfield!(model, when, res)
  elseif outboot == true
    setfield!(model, when, res[1])
    setfield!(model, :boot_estimates, res[2])
  end
  return model
end

"""

"""
function caliper!(
  model::cicmodel;
  fullmodel::Union{cicmodel, Nothing} = nothing
)

  if !isnothing(fullmodel) == true;
    model.matches = deepcopy(fullmodel.matches)
  end;

  c1 = model.stratvar == Symbol("")
  # c2 = typeof(model.caliper) .== Dict{Symbol,Float64}

  if c1
    nm, lostsi, nleft, norig = caliper(model.caliper, model.matches);
    model.matches = nm;

    model.treatedleft = nleft;
    model.treatednum = norig;

    model.matches5 = refine(nm, model.refinementnum)

    # need treatednum and treatedleft
  elseif !c1

    #= user should define each stratum caliper beforehand
    then, here, select that s from matches, and apply caliper
    to new object
    then overwrite at the end
    then refine
    =#

    nm = similar(model.matches, 0)
    sv = Symbol(string(model.stratvar) * "_stratum")

    # estblish correct types for treatedleft, treatenum
    model.treatedleft = Dict{Int64,Int64}()
    model.treatednum = Dict{Int64,Int64}()

    svals = unique(model.matches[sv])
    
    # create empties for strata without calipers, so it will plot
    for sval in svals
      if isnothing(get(model.caliper, sval, nothing))
        model.caliper[sval] = Dict{Symbol,Float64}()
      end
    end

    for (key, value) in model.caliper

      println(
        string(model.stratvar) * " stratum " * " no. " * string(key) * ":"
      )

      ms = findall(model.matches[sv] .== key);

      mcalsi, lostsi, nleft, norig = caliper(
        model.caliper[key],
        model.matches[ms, :]
      )

      append!(nm, mcalsi)

      
      model.treatedleft[key] = nleft
      model.treatednum[key] = norig
      
    end

    model.matches = nm
    model.matches5 = refine(nm, model.refinementnum)

    # @linq model.matches5 |>
    #   groupby([:tunit, :ttime]) |>
    #   combine(sum(:possible)) |>
    #   unique(:possible_function)

  else
    error("type issues, probably")
  end

  return model
end

## handle plots

function handle_attsum(model, startday, pltl)
  sumres = weeklyatt(
    model.boot_estimates, model.results, startday
  )

  lsv = (model.stratvar == Symbol(""))
  sn = pltl * namemodel(model) * "_att_wk" * ".png"
  if lsv
    pl = plot_att_sum(sumres; savename = sn)
  elseif !lsv
    pl = plot_att_sum(sumres, model.stratvar::Symbol; savename = sn)
  end
  return sumres, pl
end

function handle_attsum(model, startday, pltl, stratdict)
  sumres = weeklyatt(
    model.boot_estimates, model.results, startday
  )

  lsv = (model.stratvar == Symbol(""))
  sn = pltl * namemodel(model) * "_att_wk" * ".png"

  if !lsv
    pl = plot_att_sum(sumres, model.stratvar::Symbol, stratdict; savename = sn)
  end
  return sumres, pl
end

function handle_att!(model, pltl; post = true)
  lsv = model.stratvar == Symbol("")
  p = (post == true)

  if !p & lsv
    model.pl_att_pre = plot_att(
      model.results_pre;
      savename = pltl * namemodel(model) * "_att_pre" * ".png"
    )
  elseif !p & !lsv
    model.pl_att_pre = plot_att(model.results_pre, model.stratvar;
      savename = pltl * namemodel(model) * "_att_pre" * ".png",
      xinch = 15inch, yinch = 6inch,
      treatment = model.treatment
    )
  elseif p & lsv
    model.pl_att_post = plot_att(
      model.results;
      savename = pltl * namemodel(model) * "_att_post" * ".png"
    )
  elseif p & !lsv
    model.pl_att_post = plot_att(model.results, model.stratvar;
      savename = pltl * namemodel(model) * "_att_post" * ".png",
      xinch = 15inch, yinch = 6inch,
      treatment = model.treatment
    )
  end
  return model
end

function handle_att!(model, pltl, stratdict; post = true)
  lsv = model.stratvar == Symbol("")
  p = (post == true)

  if !p & !lsv
    model.pl_att_pre = plot_att(
      model.results_pre,
      model.stratvar,
      stratdict;
      savename = pltl * namemodel(model) * "_att_pre" * ".png",
      xinch = 15inch, yinch = 6inch,
      treatment = model.treatment
    )
  elseif p & !lsv
    model.pl_att_post = plot_att(
      model.results,
      model.stratvar,
      stratdict;
      savename = pltl * namemodel(model) * "_att_post" * ".png",
      xinch = 15inch, yinch = 6inch,
      treatment = model.treatment
    )
  end
  return model
end

function handle_balance!(model, pltl; post = true)

  lsv = model.stratvar == Symbol("")
  p = (post == true)

  if !p & lsv
    model.pl_cb_pre = plot_balance(
      model.meanbalances_pre,
      "Pre";
      savename = pltl * namemodel(model) * "_cb_pre" * ".png"
    )
  elseif !p & !lsv
    model.pl_cb_pre = plot_balance(
      model.meanbalances_pre,
      model.stratvar, "Pre";
      savename = pltl * namemodel(model) * "_cb_pre" * ".png",
      xinch = 15inch, yinch = 6inch
    )
    
  elseif p & lsv
    model.pl_cb_post = plot_balance(
      model.meanbalances_post,
      "Post";
      savename = pltl * namemodel(model) * "_cb_post" * ".png"
    )
  elseif p & !lsv
    model.pl_cb_post = plot_balance(
      model.meanbalances_post,
      model.stratvar, "Post";
      savename = pltl * namemodel(model) * "_cb_post" * ".png",
      xinch = 15inch, yinch = 6inch
    )
  end
  return model
end

function handle_balance!(model, pltl, stratdict; post = true)

  lsv = model.stratvar == Symbol("")
  p = (post == true)

  if !p & !lsv
    model.pl_cb_pre = plot_balance(
      model.meanbalances_pre,
      model.stratvar, "Pre",
      stratdict;
      savename = pltl * namemodel(model) * "_cb_pre" * ".png",
      xinch = 15inch, yinch = 6inch
    )
  elseif p & !lsv
    model.pl_cb_post = plot_balance(
      model.meanbalances_post,
      model.stratvar, "Post",
      stratdict;
      savename = pltl * namemodel(model) * "_cb_post" * ".png",
      xinch = 15inch, yinch = 6inch
    )
  end
  return model
end

"""
        runmodel(mloc)

execute script at mloc. prevents global variable scope.
"""
function runmodel(mloc::String)
    include(mloc)
end

"""
        runmodel(mlocs)

execute scripts at mlocs. prevents global variable scope.
"""
function runmodel(mloc::Vector{String})
  for ml in mloc
    include(ml)
  end
end