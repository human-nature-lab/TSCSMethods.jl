module TSCSMethods

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

funlist = [
    "types.jl",
    "construction.jl",
    "matching_setup.jl",
    "groupindices.jl",
    "retrieve_matches_utilities.jl",
    "retrieve_matches.jl",
    "retrieve_matches_missing.jl",
    "distancing_utilities.jl",
    "distancing.jl",
    "match.jl",
    "ranking.jl",
    "caliper.jl",
    "meanbalancing.jl",
    "balancing_utilities.jl",
    "overallbalancing.jl",
    "balancing.jl",
    "estimation_setup.jl",
    "estimation_utilities.jl",
    "estimation_observationweights.jl",
    "estimation.jl",
    "estimation_stratified.jl",
    "overall.jl",
    "bayesfactor.jl",
    "resampling.jl",
    "bootstrapping.jl",
    "stratification.jl",
    "refine.jl",
    "autobalancing.jl",
    "model_utilities.jl",
    "storage.jl",
    "information.jl",
    "imputation.jl",
    "inspection.jl",
    "filterunits!.jl",
    "example_data.jl"
];

for file in funlist
    include(file)
end

export
    # types
    VeryAbstractCICModel, AbstractCICModel, AbstractCICModelStratified,
    CIC, CICStratified, CaliperCIC, CaliperCICStratified,
    RefinedCIC, RefinedCICStratified, RefinedCaliperCIC, RefinedCaliperCICStratified,
    # mechanics
    match!,
    balance!, checkbalances, autobalance,
    estimate!,
    stratify, variablestrat, combostrat, customstrat,
    makemodel,
    caliper, refine,
    # saving
    name_model,
    makerecords,
    default_treatmentcategories,
    showmatches, matchinfo, obsinfo,
    save,
    # utilities
    matchprocess, quick_att, variable_filter, treatedinfo,
    relabel!, trim_model,
    # inspection utilities
    inspect_results, create_inspection_dashboard, quick_inspection, compare_strata,
    save_inspection_plots, export_results_csv,
    prepare_treatment_effect_data, prepare_counterfactual_comparison_data, prepare_percentage_change_data,
    plot_treatment_effects, plot_counterfactual_comparison, plot_percentage_changes,
    calculate_overall_summary,
    # vignette
    example_data, generate_realistic_tscs,
    filterunits!
end
