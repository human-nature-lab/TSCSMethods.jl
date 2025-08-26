# bayes–factor-t-stat.jl

# Bayesian factor calculation using R (requires RCall.jl and BayesFactor R package)
# Currently disabled to avoid R dependency in CI/CD
function bayesfactor_tstat(tstat, n)
    @warn "Bayesian factor calculation requires RCall.jl and R with BayesFactor package. Returning NaN."
    return NaN
    # Uncomment below and add RCall to dependencies to enable:
    # reval("require(BayesFactor)")
    # out = rcall(Symbol("ttest.tstat"), robject(tstat), robject(n))
    # return out[:bf][1]
end

function bfactor(b1, n)
    tstat = mean(b1) / std(b1)
    return exp(bayesfactor_tstat(tstat, n))
end
