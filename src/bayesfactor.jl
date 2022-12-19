# bayes–factor-t-stat.jl

R"""
library("BayesFactor")
"""

function bayesfactor_tstat(tstat, n)
    @rput tstat n
    
    R"""
    result <- ttest.tstat(
        t = tstat,
        n,
        nullInterval = NULL,
        rscale = "medium",
        complement = FALSE,
        simple = FALSE
    )
    """

    R"""
    lbf1 <- result[['bf']]
    """

    return @rget lbf1
end

function bfactor(b1, n)
    tstat = mean(b1) / std(b1)
    return exp(bayesfactor_tstat(tstat, n))
end
