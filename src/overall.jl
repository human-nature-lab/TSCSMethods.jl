# overall.jl

mutable struct Overall
    att::Float64
    percentiles::Vector{Float64}
    bayesfactor::Float64
    ntreatedmean::Float64
    pvalue::Float64
end

function overall(
    ; att = NaN, percentiles = Float64[], bayesfactor = NaN,
    ntreatedmean = 0.0, pvalue = NaN
)

    return Overall(
        att, percentiles, bayesfactor, ntreatedmean, pvalue
    )
end
