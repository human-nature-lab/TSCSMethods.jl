# epi_model.jl

using tscsmethods, COVIDPoliticalEvents
using DataFrames, DataFramesMeta, Dates
using JLD2:save_object,load_object

spth = model_path = ""; # "primary-elections/";

dat = load_object("/Users/emf/Library/Mobile Documents/com~apple~CloudDocs/Yale/yale research/COVID19/covid-19-data/data/cvd_dat.jld2");

# base model

cc = deathmodel("primary epi", :primary, :epi);

iter = 10000;
refinementnum = 5;

dat = dataprep(
  dat, cc;
  t_start = 0
);

# @time match!(cc, dat; distances = true);

save_object(name_model(cc), cc)

save_modelset(
  spth * name_model(cc),
  cc
);
  
cc, ccr, cal, calr, labels = load_modelset(spth * name_model(cc) * ".jld2")

#=
@time balance!(cc, dat);

@time estimate!(cc, dat; iter = iter);

ccr = make_refined(cc; refinementnum = refinementnum);

@time estimate!(ccr, dat; iter = iter);

## caliper model

caliper = Dict(
  covariates[1] => 0.25,
  covariates[2] => 0.25,
  covariates[3] => 0.25
);

@time cal = make_caliper(cc, caliper);

@time estimate!(cal, dat; iter = iter);

@time calr = make_refined(cal; refinementnum = refinementnum);

@time estimate!(calr, dat; iter = iter);

# name and save models
save_modelset(
  spth * name_model(cc),
  cc; ccr = ccr, cal = cal, calr = calr
  # labels = labels
);

plot_modelset(
  model_path;
  variablecolors = varcol, # from paper pkg.
  base_savepath = "" # ends in / if not blank
);
=#