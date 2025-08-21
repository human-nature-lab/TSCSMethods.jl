# dependencies.jl

using LinearAlgebra:permutedims, pinv, diag, diagind
using DataFrames, DataFramesMeta
using StatsBase:cov,mean,std,sample,sample!
using Statistics:quantile
using Distances:mahalanobis, weuclidean
using Parameters
using Accessors:@set,@reset
using FLoops
using Random
using Dates
import JLD2:save_object,load_object
import CSV

# Optional R dependency for Bayesian factors (disabled by default)
# using RCall
# R dependencies: BayesFactor
