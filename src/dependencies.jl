# dependencies.jl

using LinearAlgebra:permutedims, pinv, diag, diagind
using DataFrames, DataFramesMeta
using StatsBase:cov,mean,std,sample
using Statistics:quantile
using Distances:mahalanobis, weuclidean
using Parameters
using Accessors:@set,@reset
using FLoops

using CairoMakie
using Colors
