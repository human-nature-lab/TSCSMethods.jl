# model.jl

"""
    variable_filter(
      model, variable, dat;
      mn = nothing, mx = nothing
    )

Remove treated observations according to some variable values.
"""
function variable_filter(
  model, variable, dat;
  mn = nothing, mx = nothing
)

  @unpack t, id, treatment = model

  # remove elections prior to March 10

  dt = @subset(dat, $treatment .== 1);

  tple = [(dt[i, t], dt[i, id]) for i in eachindex(dt[!, t])];
  dict = Dict(tple .=> dt[!, variable]);

  obinclude = fill(false, length(model.observations));
  for (i, ob) in enumerate(model.observations)
    cond = true
    if !isnothing(mn)
      cond = cond & (dict[ob] >= mn)
    end
    if !isnothing(mx)
      cond = cond & (dict[ob] <= mx)
    end
      obinclude[i] = cond
  end

  @reset model.observations = model.observations[obinclude];
  @reset model.matches = model.matches[obinclude];
  # @reset model.results = TSCSMethods.DataFrame();

  @reset model.treatednum = length(model.observations)
  
  return model
end

"""
    function treatedinfo(
      model, variables, dat;
    )

Gives variable values for treated observations present in the model, for
the chosen set of variables. Order is the same as model.observations.
"""
function treatedinfo(
  model, variables, dat;
)

  @unpack t, id, treatment = model

  # remove elections prior to March 10

  dt = @subset(dat, $treatment .== 1);

  tple = [(dt[i, t], dt[i, id]) for i in eachindex(dt[!, t])];
  
  obout = DataFrame(
    obs = model.observations,    
    )
    
    for variable in variables
      obout[!, variable] = Vector{eltype(dat[!, variable])}(undef, nrow(obout))
      dict = Dict(tple .=> dt[!, variable]);
      
      for (i, e) in enumerate(obout.obs)
        obout[i, variable] = dict[e]
      end
    end
  return obout
end

function relabel!(
  calmodel, refcalmodel, dat; stratifier = nothing, digits = 2
)

  stratifier = if isnothing(stratifier)
    calmodel.stratifier
  else
    stratifier
  end

  # the cal model and refcalmodel have the same stratum ranges
  calinfo = treatedinfo(
    calmodel, [stratifier], dat;
  )
  calinfo[!, :stratum] = calmodel.strata

  calinfo = @chain calinfo begin
    groupby(:stratum)
    combine(stratifier => extrema => stratifier)
  end

  exts = Dict(calinfo.stratum .=> calinfo[!, stratifier])

  relabels = Dict{Int, String}();
  for s in sort(collect(keys(exts)))
    mn, mx = exts[s]
    if typeof(mn) == Float64
      mn = round(mn; digits = digits)
      mx = round(mx; digits = digits)
    end

    relabels[s] = if ismissing(mn) & ismissing(mx)
      "Missing Values"
    else
      if mn < mx
        string(mn) * " to " * string(mx)
      elseif mn == mx
        string(mn)
      end
    end
  end

  for (k,v) in relabels; calmodel.labels[k] = v end # add labels
  for (k,v) in relabels; refcalmodel.labels[k] = v end # add labels
  return calmodel, refcalmodel
end

function relabel!(
  m, dat; stratifier = nothing, digits = 2
)

  stratifier = if isnothing(stratifier)
    m.stratifier
  else
    stratifier
  end

  # the cal model and refcalmodel have the same stratum ranges
  calinfo = treatedinfo(
    m, [stratifier], dat;
  )
  calinfo[!, :stratum] = m.strata

  calinfo = @chain calinfo begin
    groupby(:stratum)
    combine(stratifier => extrema => stratifier)
  end

  exts = Dict(calinfo.stratum .=> calinfo[!, stratifier])

  relabels = Dict{Int, String}();
  for s in sort(collect(keys(exts)))
    mn, mx = exts[s]
    if typeof(mn) == Float64
      mn = round(mn; digits = digits)
      mx = round(mx; digits = digits)
    end
    relabels[s] = if ismissing(mn) & ismissing(mx)
      "Missing Values"
    else
      if mn < mx
        string(mn) * " to " * string(mx)
      elseif mn == mx
        string(mn)
      end
    end
  end

  for (k,v) in relabels; m.labels[k] = v end # add labels
  return m
end
