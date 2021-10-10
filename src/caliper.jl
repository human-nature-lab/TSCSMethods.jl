# caliper.jl

"""
this should be written to simply add the exclude column

ranking / counting should be a separate function call
"""
function applycaliper(cc::AbstractCICModel, covcaliper)

  chk = Matrix{Bool}(undef, nrow(cc.matches), length(cc.covariates));
  for (c, covar) in enumerate(cc.covariates)
      chk[:, c] = cc.matches[!, covar] .<= covcaliper[covar]
  end

  return [all(r) for r in eachrow(chk)]
end

#=
  for covar in cc.covariates
    cc.matches[!, covar] .>= cc.caliper[covar]
  end

  cc.matches[!, :exclude] .= false;
  idx = Int64[];
  for covar in cc.covariates
    append!(idx, findall(cc.matches[!, covar] .>= cc.caliper[covar]))
  end

  idx = sort(unique(idx));
  cc.matches[idx, :exclude] .= true

  cc_cal.matches = @subset(cc.matches, :exclude .== false)
  rank!(cc_cal::cicmodel)

  return cc_cal
end
=#

function rank!(cc::AbstractCICModel)
  cc.matches = @chain cal.matches begin
    @orderby(:treattime, :treatunit, :f, :mdist)
    groupby([:treattime, :f, :treatunit])
    transform( 
      :matchunit => eachindex => :rank,
      nrow => :numpossible_cal # name should depend on whether it is caliper type thing or not
    )
  end
  return cc
end

#= balance caliper

#=
figure out which matches to drop based on balance score
1. check over the whole fmin + mmin to fmax + mmax range
2. record position in length 80 vector that must be dropped
based on the position, we know which f(s) to remove from matches

this is very strict, stating that no match may have a std. difference above cutoff at any point
=# 

# generalize cutoff to a dict for each covariates
function flag_balances(
  cc::cicmodel,
  cutoff::Float64 = 1/4
)

  dropwhere = Dict{Tuple{Int64, Int64, Int64}, Vector{Int64}}()

  begin
    for i in 1:nrow(cc.balances)
      x = Vector{Int64}()
      for covar in cc.covariates
        if !cc.timevary[covar]
          if (abs(@views cc.balances[i, covar]) > cutoff)
          
            dropwhere[
                cc.balances[i, :treattime],
                cc.balances[i, :treatunit],
                cc.balances[i, :matchunit],
              ] = collect(1:80)
            continue
          end
        else
          append!(x, findall(
            skipmissing(
              (abs.(@views cc.balances[i, covar]) .> cutoff)
            )
          ));
        end
        if length(x) > 0
          dropwhere[
              cc.balances[i, :treattime],
              cc.balances[i, :treatunit],
              cc.balances[i, :matchunit],
            ] = x
        end
      end
    end
  end
  return dropwhere
end

function flag_balances(cc::cicmodel)

  dropwhere = Dict{Tuple{Int64, Int64, Int64}, Union{Vector{Int64}, Bool}}()

  begin
    for i in 1:nrow(cc.balances)
      x = Vector{Int64}()
      for covar in cc.covariates
        if !cc.timevary[covar]
          if (abs(@views cc.balances[i, covar]) > cc.caliper[covar])
          
            dropwhere[
                cc.balances[i, :treattime],
                cc.balances[i, :treatunit],
                cc.balances[i, :matchunit],
              ] = true
            continue
          end
        else
          append!(x, findall(
            skipmissing(
              (abs.(@views cc.balances[i, covar]) .> cc.caliper[covar])
            )
          ));
        end
        if length(x) > 0
          dropwhere[
              cc.balances[i, :treattime],
              cc.balances[i, :treatunit],
              cc.balances[i, :matchunit],
            ] = x
        end
      end
    end
  end
  return dropwhere
end

"""
use dropwhere to flag matches
inputs dropwhere, from flag balances

adds :exclude col. to matches
"""
function flag_matches!(cc::cicmodel, dropwhere)
  lmm = length(cc.mmin + cc.fmin:cc.mmax + cc.fmax)
  @eachrow! cc.matches begin
    @newcol :exclude::Vector{Bool}
    post = get(dropwhere, (:treattime, :treatunit, :matchunit), Vector{Int64}())
    if length(post) == lmm
      :exclude = true
    elseif length(post) > 0
      if any((post .> (cc.mmin + :f)) .& (post .< (cc.mmax + :f))) # any point in matching over thresh excludes that match
        :exclude = true
      end
    else :exclude = false
    end
  end
  return cc
end

"""
apply match-level caliper to standardized balance scores

"""
function balance_caliper!(cc::cicmodel)

  dropwhere = flag_balances(cc)
  flag_matches!(cc, dropwhere)

  cc.matches = cc.matches[.!cc.matches[!, :exclude], :];

  # keep the balances as they are average balance calculation will exclude
  # cc.balances = leftjoin(
  #   unique(cc.matches[!, [:treattime, :treatunit, :matchunit]]),
  #   cc.balances,
  #   on = [:treattime, :treatunit, :matchunit]
  # );

  rank!(cc)

  return cc
end
=#