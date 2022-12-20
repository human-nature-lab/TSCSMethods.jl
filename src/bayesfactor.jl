# bayes–factor-t-stat.jl

function bayesfactor_tstat(tstat, n)
    reval("require(BayesFactor)")
    
    out = rcall(Symbol("ttest.tstat"), robject(tstat), robject(n))
    return out[:bf][1]
end

function bfactor(b1, n)
    tstat = mean(b1) / std(b1)
    return exp(bayesfactor_tstat(tstat, n))
end
