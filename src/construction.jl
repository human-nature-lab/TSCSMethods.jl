# construction.jl

"""
    make_caliper(cc::cicmodel, caliper)

Create the caliper model from the full cic model, calculates the average, and overall balance scores, gives the number of treated observations that survive the caliper.
"""
function make_caliper(cc::cicmodel, caliper)
  # caliper::Dict{Symbol, Float64}()
  inc = applycaliper(cc, caliper);

  cal = calipercicmodel(
    title = cc.title * " caliper",
    id = cc.id,
    t = cc.t,
    outcome = cc.outcome,
    treatment = cc.treatment,
    covariates = cc.covariates,
    timevary = cc.timevary,
    fmin = cc.fmin,
    fmax = cc.fmax,
    stratifier = cc.stratifier,
    reference = cc.reference,
    mmin = cc.mmin,
    mmax = cc.mmax,
    matches = @view(cc.matches[inc, :]), # caliper here
    balances = @view(cc.balances[!, :]),
    # meanbalances = , # do func
    # grandbalances = , # do func
    iterations = cc.iterations,
    results = DataFrame(),
    treatednum = Int64(0),
    treatedleft = Int64(0),
    estimator = "ATT",
    caliper = caliper,
    fullmod = Base.RefValue{cicmodel}
  )

  if all(.!inc)
    return "No matches survive the caliper(s)."
  end

  meanbalance!(cal);

  # meanbalance! does not carry over the stratum label, so reassign
  if cc.stratifier != Symbol("") 
    uset = unique(
      cc.meanbalances[!, [:treattime, :treatunit, :stratum]],
      view = true
    )
    cal.meanbalances = leftjoin(
      cal.meanbalances, uset, on = [:treattime, :treatunit]
    )
  end

  grandbalance!(cal);

  cal.treatednum = cc.treatednum;
  if cc.stratifier == Symbol("")
    cal.treatedleft = nrow(unique(cal.meanbalances[!, [:treattime, :treatunit]]));
  else
    treatedlefts!(cal) # left
  end

  trtxt = string(cal.treatedleft) * " treated observations remain.";
  tntxt = "(" * "out of " * string(cc.treatednum) * ")";
  println(trtxt * tntxt)

  return cal
end

"""
    make_refined(cc::AbstractCICModel; refinementnum = 5)

Refine a full or caliper model less than or equal to the chosen number of best matches per treated observation. Calculates the average, and overall balance scores, gives the number of treated observations that survive the caliper.
"""
function make_refined(cc::AbstractCICModel; refinementnum = 5)

  if typeof(cc) == calipercicmodel
    calip = cc.caliper
  else
    calip = Dict{Symbol, Float64}()
  end;

  rf = refinedcicmodel(
    title = cc.title * " refined",
    id = cc.id,
    t = cc.t,
    outcome = cc.outcome,
    treatment = cc.treatment,
    covariates = cc.covariates,
    timevary = cc.timevary,
    fmin = cc.fmin,
    fmax = cc.fmax,
    stratifier = cc.stratifier,
    reference = cc.reference,
    mmin = cc.mmin,
    mmax = cc.mmax,
    matches = refine(cc, refinementnum),
    balances = @view(cc.balances[!, :]),
    # meanbalances = , # do func
    # grandbalances = , # do func
    iterations = cc.iterations,
    results = DataFrame(),
    treatednum = Int64(0),
    treatedleft = Int64(0),
    estimator = "ATT",
    caliper = calip,
    fullmod = Base.RefValue{cicmodel}
  )

  if nrow(cc.matches) == 0
    return "There are no matches."
  end

  meanbalance!(rf)
  
  # meanbalance! does not carry over the stratum label, so reassign
  if cc.stratifier != Symbol("") 
    uset = unique(
      cc.meanbalances[!, [:treattime, :treatunit, :stratum]],
      view = true
    )
    rf.meanbalances = leftjoin(
      rf.meanbalances, uset, on = [:treattime, :treatunit]
    )
  end
  
  grandbalance!(rf)

  rf.treatednum = cc.treatednum;
  if cc.stratifier == Symbol("")
    rf.treatedleft = nrow(unique(cc.meanbalances[!, [:treattime, :treatunit]]));
  else
    treatedlefts!(rf)
  end

  return rf
end
