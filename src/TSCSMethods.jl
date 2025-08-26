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
    "core/types.jl",
    "core/construction.jl",
    "matching/matching_setup.jl",
    "advanced/groupindices.jl",
    "matching/retrieve_matches_utilities.jl",
    "matching/retrieve_matches.jl",
    "matching/retrieve_matches_missing.jl",
    "matching/distancing_utilities.jl",
    "matching/distancing.jl",
    "matching/match.jl",
    "matching/ranking.jl",
    "matching/caliper.jl",
    "balancing/meanbalancing.jl",
    "balancing/balancing_utilities.jl",
    "balancing/overallbalancing.jl",
    "balancing/balancing.jl",
    "estimation/estimation_setup.jl",
    "estimation/estimation_utilities.jl",
    "estimation/estimation_observationweights.jl",
    "estimation/estimation.jl",
    "estimation/estimation_stratified.jl",
    "estimation/overall.jl",
    "estimation/bayesfactor.jl",
    "estimation/resampling.jl",
    "estimation/bootstrapping.jl",
    "advanced/stratification.jl",
    "advanced/refine.jl",
    "balancing/autobalancing.jl",
    "utilities/model_utilities.jl",
    "utilities/storage.jl",
    "utilities/information.jl",
    "utilities/imputation.jl",
    "utilities/inspection.jl",
    "utilities/filterunits!.jl",
    "utilities/example_data.jl"
];

for file in funlist
    include(file)
end

export
    # types
    VeryAbstractCICModel, AbstractCICModel, AbstractCICModelStratified,
    CIC, CICStratified, CaliperCIC, CaliperCICStratified,
    RefinedCIC, RefinedCICStratified, RefinedCaliperCIC, RefinedCaliperCICStratified,
    TreatmentObservationMatches, Overall,
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
