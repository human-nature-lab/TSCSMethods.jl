# overall.jl

mutable struct Overall
    att::Union{Float64, Vector{Float64}}
    percentiles::Union{Vector{Float64}, Vector{Vector{Float64}}}
    bayesfactor::Union{Float64, Vector{Float64}}
    ntreatedmean::Union{Float64, Vector{Float64}}
    pvalue::Union{Float64, Vector{Float64}}
    stratum::Union{Missing, Vector{Float64}}
end

function overall(
    ;
    att = NaN,
    percentiles = Float64[],
    bayesfactor = NaN,
    ntreatedmean = 0.0,
    pvalue = NaN,
    stratum = missing
)

    return Overall(
        att, percentiles, bayesfactor, ntreatedmean, pvalue, stratum
    )
end
