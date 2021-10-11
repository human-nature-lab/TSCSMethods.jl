# example
# without package format, loads functions directly

using Random
Random.seed!(2019)

# cvd_polev only
import CSV
using Dates

# dev
using JLD2:load_object
using LinearAlgebra:permutedims, pinv, diag, diagind
using DataFrames, DataFramesMeta
using StatsBase:cov,mean,std,sample
using Statistics:quantile
using Distances:mahalanobis, weuclidean
using Parameters:@with_kw

dat = load_object("/Users/emf/Library/Mobile Documents/com~apple~CloudDocs/Yale/yale research/COVID19/covid-19-data/data/cvd_dat.jld2");

id = :fips; t = :running;
mmin = -50; mmax = -1;
fmin = 10; fmax = 40;
treatment = :primary; outcome = :death_rte;
covariates = calvars = [Symbol("Pop. Density"), :death_rte];

@subset!(dat, :running .> 0, :date .< Date("2020-10-01"));

covariates = [
  Symbol("Pop. Density"),
  Symbol("Cum. Death Rate"),
  Symbol("Date of First Case")
];

mmin = -50; mmax = -1; tpoint = -1;

# problem here
# using tscsmethods:cicmodel

include("pkg_types.jl")
include("matching.jl")
include("balancing.jl")
include("caliper.jl")
include("estimation.jl")

cc = load_object("cc.jld2");

#=
cc = cicmodel(
  title = "test model",
  id = id,
  t = t,
  outcome = outcome,
  treatment = treatment,
  covariates = covariates,
  timevary = Dict(
    covariates[1] => false, covariates[2] => true, covariates[3] => false
  ),
  fmin = fmin,
  fmax = fmax,
  mmin = mmin,
  mmax = mmax,
  matches = cc.matches,
  balances = cc.balances,
  meanbalances = cc.meanbalances
);

cc.matches = deepcopy(full.matches);
cc.balances = deepcopy(full.balances);

# import JLD2; JLD2.save_object("cc.jld2", cc);
=#

# match

# @time match!(cc, dat; distances = true);

# balance

# @time balance!(cc, dat);

@time chk = balancecheck(cc);

# refine
ccr = make_refined(cc; refinementnum = 5);

# broken in not-strat case
# ccr_chk = balancecheck(ccr);

# caliper

caliper = Dict(
  covariates[1] => 0.5,
  covariates[2] => 0.5,
  covariates[3] => 0.5
);

@time cal = make_caliper(cc, caliper);

# refine caliper

@time calr = make_refined(cal; refinementnum = 5);

# estimation

@time estimate!(cc, dat; iter = 500);

@time estimate!(ccr, dat; iter = 500);

@time estimate!(cal, dat; iter = 500);

@time estimate!(calr, dat; iter = 500);

# stratification example
# (start with balanced and matched model cicmodel)

# include("stratified.jl");
include("stratification.jl");

@time cc, labels = stratify!(datestrat!, [cc]);

# refine
ccr = make_refined(cc; refinementnum = 5);

ccr_chk = balancecheck(ccr);

# caliper

caliper = Dict(
  covariates[1] => 0.5,
  covariates[2] => 0.5,
  covariates[3] => 0.5
);

@time cal = make_caliper(cc, caliper);

# refine caliper

@time calr = make_refined(cal; refinementnum = 5);

# estimation

@time estimate!(cc, dat; iter = 500);

@time estimate!(ccr, dat; iter = 500);

@time estimate!(cal, dat; iter = 500);

@time estimate!(calr, dat; iter = 500);

# plotting

using Colors, CairoMakie

# while modeling
plot_cb(
  cc;
  labels = labels,
  variablecolors = varcol,
  fw = 700, fl = 300,
  spath = nothing
)

plot_cbs(
  cc,
  ccr;
  labels = labels,
  variablecolors = varcol,
  fw = 700,
  fl = 300,
  spath = nothing
)
#

# paper plot
mp = model_pl(
  cc;
  labels = labels
);

# name and save models
save_modelset(
  spth * name_model(cc::AbstractCICModel),
  cc; ccr = ccr, cal = cal, calr = calr,
  labels = labels
);

plot_modelset(
  model_path;
  variablecolors = varcol, # from paper pkg.
  base_savepath = "" # ends in /
);
