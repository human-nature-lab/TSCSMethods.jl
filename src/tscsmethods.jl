module tscsmethods

# load all dependencies

using StatsBase, Statistics

import LinearAlgebra.Diagonal
import IterTools.product

using Gadfly
import Cairo, Fontconfig

using CSV, DataFrames, DataFramesMeta

# include code files here
include("helper_functions.jl")
include("matching_functions.jl")
include("balancing_functions.jl")
include("estimation_functions.jl")
include("get_wits.jl")
include("bootstrap_estimation.jl")
include("plotting_functions.jl")
# include("load_and_clean.jl"); doesn't go into pkg
# need to explain how the data should be formatted, etc.


end # module tscsmethods
