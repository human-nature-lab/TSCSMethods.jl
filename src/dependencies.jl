# dependencies.jl

using CSV, DataFrames, DataFramesMeta
using StatsBase, Statistics, Gadfly, Compose, Parameters

import LinearAlgebra.Diagonal, LinearAlgebra

import IterTools.product

import Cairo, Fontconfig

import JLD2, Dates