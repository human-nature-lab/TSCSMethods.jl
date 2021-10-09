# dependencies.jl

using JLD2:load_object,save_object
using LinearAlgebra:permutedims, pinv, diag, diagind

using DataFrames, DataFramesMeta
using StatsBase:cov,mean,std,sample
using Statistics:quantile
using Distances:mahalanobis, weuclidean
using Parameters:@with_kw

using CairoMakie
import Colors
