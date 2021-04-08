# model_functions.jl
#= 
  set of functions needed to execute specific types of models, observe
  model properties, etc.
=#

import JLD2.@save, JLD2.@load

"""
      datamake(covariates, varset, fmin, fmax, tmin,
      id, tpoint, outcome, treatment
      )

use this function when the original data is mutated in some model run...
"""
function datamake(
  covariates, varset,
  fmin, fmax, tmin,
  id, tpoint, outcome, treatment
)

  @load "../covid-19-data/data/cvd_dat.jld2"

  dat = @where(dat, :running .>= 0);

  dat = remove_incomplete(
    dat, covariates, 9, 40, tmin, id, t, outcome, tpoint);

  select!(dat, varset);
  return dat
end

function evthres!(
    dat,
    ct::Int64, ot::Int64,
    fmin::Int64, fmax::Int64,
    sz::Symbol, trt::Symbol, t::Symbol, id::Symbol
  )
  
    # size <= ct -> control (switch treated to zero)
    # size > ct & size <= ot -> omit (remove treated day to 40 days after from data)
    # N.B. we will now have to handle dropped observations in the data
    # - CHECK THIS
  
    # this occurs prior to matching, and changes require new matching
    # output a new dat object (don't overwrite) -- use referencing in future vers.
  
    # probably a smarter way to do this via Ref / @views / etc:
    # actually, just mutate
    # dat = deepcopy(dat);
  
    rning = dat[!, t];
    idv = dat[!, id];
    trv = dat[!, trt];
    sv = dat[!, sz];
  
    # create new treatment variable
    # trnnme = Symbol(String(trt) * "thr");
    # trn = copy(trv);
  
    nix = Vector{Int64}(undef, 0);
  
    for i in 1:length(trv)
      trvi = @views(trv[i]);
      idi = @views(idv[i]);
      ti = @views(rning[i]);
      svi = @views(sv[i]);
      if trvi == 1
        if svi <= ct
          trv[i] = 0
        elseif (svi > ct) & (svi <= ot)
          # hard part
          # remove data for that unit from treatment to treatment + fmax
          # easiest way: use findall, also safest
          anix = findall(
            (idv .== idi) .& (rning .>= ti) .& (rning .<= (ti + fmax))
          )
          append!(nix, anix)
        end
      end
    end
    deleterows!(dat, sort(unique(nix)))
    return dat
  end

"""
remove some set of units, meant for adj removal -- omission when worried about
spillover to some set of counties
"""
function removeadj!(dat::DataFrame, remset::DataFrame; offset = 40)
    RDX = Vector{Int64}();
    for i = eachindex(remset.fips)
        c1 = dat[!, id] .== remset[id][i];
        c2 = dat[!, t] .>= remset[t][i];
        c3 = dat[!, t] .<= (remset[t][i] + offset);
        append!(RDX, findall(c1 .& c2 .& c3))
    end
    return delete!(dat, unique(sort(RDX)))
end

"""
get the properties for the set of treated units in a model
"""
function gettreatedproperties(model::cicmodel, t, id, dat::DataFrame)
    tobs = @where(model.matches5, :tunit .== :munit)
  
    idx = Vector{Int64}();
    for i = eachindex(tobs.ttime)
      append!(
        idx,
        getvarind(dat, @views(tobs.tunit[i]), id, t, @views(tobs.ttime[i]))
      )
    end
    return dat[idx, :]
  end
  