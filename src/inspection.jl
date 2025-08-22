# inspection.jl
#
# TSCSMethods Inspection and Visualization Utilities
#
# This module provides comprehensive tools for inspecting and visualizing results from 
# TSCSMethods causal inference analyses. The architecture separates data processing from 
# visualization to enable flexible, reusable, and testable code.
#
# ## Quick Start Example
# ```julia
# using TSCSMethods
# 
# # Run TSCSMethods analysis
# data = example_data()
# model = makemodel(data, :t, :id, :gub, :Y, [:X1, :X2], 
#                   Dict(:X1 => false, :X2 => false), -15:-10, 1:5)
# match!(model, wids=50)
# balance!(model)  
# estimate!(model, dobayesfactor=false)
#
# # Generate counterfactual analysis
# imputation = impute_results(model, matches, data, :t, :id)
#
# # Comprehensive inspection
# inspection = inspect_results(imputation, :Y)
# inspection.treatment_effects_plot     # ATT over time
# inspection.counterfactual_plot        # Observed vs counterfactual
# inspection.summary                    # Summary statistics
#
# # Or create dashboard
# dashboard = create_inspection_dashboard(imputation, :Y)
# ```
#
# ## Architecture Overview
#
# ### 1. Data Processing Functions (Pure, No Plotting Dependencies)
# - `prepare_treatment_effect_data()` - Extract ATT data with confidence intervals
# - `prepare_counterfactual_comparison_data()` - Format observed vs counterfactual data
# - `prepare_percentage_change_data()` - Extract percentage change data
# - `calculate_overall_summary()` - Compute summary statistics
#
# ### 2. Plotting Functions (Take Prepared Data, Handle Visualization Only)  
# - `plot_treatment_effects()` - ATT over time with confidence intervals
# - `plot_counterfactual_comparison()` - Trajectory comparison with confidence bands
# - `plot_percentage_changes()` - Treatment effects as percentage changes
#
# ### 3. High-Level Inspection Functions (Coordinate Everything)
# - `inspect_results()` - Main inspection function with clean API
# - `create_inspection_dashboard()` - Comprehensive multi-panel dashboard
#
# ## Dependencies
# - **Always available**: Data processing functions work with base Julia + DataFrames
# - **Optional**: Plotting functions require CairoMakie.jl (install with `Pkg.add("CairoMakie")`)
# - **Graceful fallback**: Clear error messages when plotting dependencies unavailable

# Import required functions for data processing (already available from main module)
using Statistics: mean
using DataFrames: names
using Base: merge  # For merging NamedTuples

# Check for plotting dependencies  
const PLOTTING_AVAILABLE = try
    @eval using CairoMakie
    true
catch
    @warn "Plotting functionality requires CairoMakie.jl. Install with: Pkg.add(\"CairoMakie\")"
    false
end

## Data Processing Functions (separated from plotting)

"""
    prepare_treatment_effect_data(imputation_results::ImputationResults)

Extract and sort treatment effect data for plotting.

This function processes ImputationResults to create a clean DataFrame containing
treatment effects (ATT) and confidence intervals organized by time period.

# Arguments
- `imputation_results::ImputationResults`: Output from `impute_results()`

# Returns
- `DataFrame`: Sorted by time period, containing columns:
  - `f`: Time periods  
  - `att`: Average treatment effects for each time period
  - `2.5%`: Lower bounds of 95% confidence intervals  
  - `97.5%`: Upper bounds of 95% confidence intervals

# Example
```julia
imputation = impute_results(model, matches, data, :t, :id)
results_df = prepare_treatment_effect_data(imputation)

# Access the data (now just DataFrame columns)
plot_times = results_df.f           # [-4, -3, -2, -1, 1, 2, 3, 4, 5]
effects = results_df.att            # [0.1, 0.2, 0.15, 0.05, 0.3, 0.4, ...]
lower_ci = results_df[!, Symbol("2.5%")]  # [-0.1, 0.0, -0.05, -0.15, 0.1, ...]
```
"""
function prepare_treatment_effect_data(imputation_results::ImputationResults)
    # Return sorted DataFrame - much simpler!
    return sort(imputation_results.results, :f)
end

"""
    prepare_counterfactual_comparison_data(imputation_results::ImputationResults, outcome_var::Symbol)

Extract and sort data for observed vs counterfactual trajectory plotting.

This function prepares the DataFrame needed to visualize the key causal inference comparison:
what actually happened (observed) vs. what would have happened without treatment (counterfactual).

# Arguments
- `imputation_results::ImputationResults`: Output from `impute_results()`
- `outcome_var::Symbol`: Column name of the outcome variable (e.g., `:Y`, `:deaths`, `:GDP`)

# Returns
- `DataFrame`: Sorted by time period, containing all original columns plus outcome variable.
  Key columns for plotting:
  - `f`: Time periods
  - `[outcome_var]`: Actual observed outcomes for treated units  
  - `counterfactual_trajectory`: Estimated counterfactual outcomes
  - `counterfactual_lower`: Lower bounds of counterfactual confidence intervals
  - `counterfactual_upper`: Upper bounds of counterfactual confidence intervals
  - `att`: Treatment effects (observed - counterfactual)

# Example
```julia
imputation = impute_results(model, matches, data, :t, :id)  
results_df = prepare_counterfactual_comparison_data(imputation, :Y)

# The key causal story (now just DataFrame columns)
observed_trajectory = results_df[!, :Y]                    # What actually happened
counterfactual_trajectory = results_df.counterfactual_trajectory  # What would have happened
treatment_effect = results_df.att                          # The causal impact
```

# Notes
The counterfactual trajectory represents the estimated outcome path that treated units
would have followed if they had not received treatment, constructed by averaging
outcomes from matched control units.
"""
function prepare_counterfactual_comparison_data(imputation_results::ImputationResults, outcome_var::Symbol)
    # Validate outcome variable exists
    outcome_var ∉ names(imputation_results.results) && throw(ArgumentError("Outcome variable $outcome_var not found in results"))
    
    # Return sorted DataFrame - outcome_var column should already exist
    return sort(imputation_results.results, :f)
end

"""
    prepare_percentage_change_data(imputation_results::ImputationResults)

Extract and sort percentage change data for plotting.

# Arguments  
- `imputation_results::ImputationResults`: Output from `impute_results()`

# Returns
- `DataFrame`: Sorted by time period, containing columns:
  - `f`: Time periods
  - `pct_change`: Treatment effects as percentage changes
  - `pct_change_lower`: Lower bounds of percentage change confidence intervals  
  - `pct_change_upper`: Upper bounds of percentage change confidence intervals

# Example
```julia
imputation = impute_results(model, matches, data, :t, :id)
results_df = prepare_percentage_change_data(imputation)

# Access percentage data (now just DataFrame columns)
pct_effects = results_df.pct_change        # [2.1, 4.2, 1.5, 0.5, ...]
pct_lower = results_df.pct_change_lower    # [0.1, 2.0, -0.5, -1.5, ...]
```
"""
function prepare_percentage_change_data(imputation_results::ImputationResults)
    # Return sorted DataFrame - same as the others!
    return sort(imputation_results.results, :f)
end

"""
    calculate_overall_summary(imputation_results::ImputationResults, outcome_var::Symbol)

Calculate overall treatment effect summary statistics and balance diagnostics.

This function computes key summary statistics that help assess the quality and
magnitude of the causal inference results.

# Arguments
- `imputation_results::ImputationResults`: Output from `impute_results()`
- `outcome_var::Symbol`: Column name of the outcome variable

# Returns
- `NamedTuple`: Summary statistics with fields:
  - `treated_observed_mean`: Average observed outcome for treated units
  - `overall_att`: Overall average treatment effect
  - `counterfactual_mean`: Average counterfactual outcome 
  - `baseline_difference`: Pre-treatment difference between treated and matched units
  - `treated_pretreatment_avg`: Pre-treatment average for treated units
  - `matched_pretreatment_avg`: Pre-treatment average for matched controls

# Example
```julia
imputation = impute_results(model, matches, data, :t, :id)
summary = calculate_overall_summary(imputation, :Y)

# Key statistics
println("Overall ATT: ", summary.overall_att)
println("Baseline balance: ", summary.baseline_difference)  # Should be close to 0
println("Treatment magnitude: ", summary.overall_att / summary.treated_observed_mean * 100, "%")
```

# Notes
- `baseline_difference` close to 0 indicates good matching quality
- `overall_att` is the main causal estimate  
- Comparing pre-treatment averages helps assess whether matching achieved balance
"""
function calculate_overall_summary(imputation_results::ImputationResults, outcome_var::Symbol)
    results_df = imputation_results.results
    
    # Overall averages
    treated_avg = mean(results_df[!, outcome_var])
    overall_att = mean(results_df.att)
    
    # Counterfactual average
    counterfactual_avg = treated_avg - overall_att
    
    summary = (
        treated_observed_mean = treated_avg,
        overall_att = overall_att,
        counterfactual_mean = counterfactual_avg,
        baseline_difference = imputation_results.baseline_difference,
        treated_pretreatment_avg = imputation_results.treated_pretreatment_avg,
        matched_pretreatment_avg = imputation_results.matched_pretreatment_avg
    )
    
    return summary
end

## High-Level Inspection Functions

"""
    inspect_results(imputation_results::ImputationResults, outcome_var::Symbol; 
                   plot_percentage_changes::Bool = false)

Comprehensive inspection of TSCSMethods results with clean separation of data processing and plotting.

# Arguments
- `imputation_results`: ImputationResults struct from impute_results()
- `outcome_var`: Symbol representing the outcome variable column name
- `plot_percentage_changes`: Whether to include percentage change plots (default: false)

# Returns
- `NamedTuple` with fields:
  - `treatment_effects_plot`: Figure showing ATT over time with confidence intervals
  - `counterfactual_plot`: Figure showing observed vs counterfactual trajectories  
  - `summary`: NamedTuple with summary statistics and balance diagnostics
  - `data`: Prepared data structures for custom plotting or further analysis
  - `percentage_plot`: (optional) Figure showing treatment effects as percentage changes

# Example
```julia
using TSCSMethods

# Complete workflow
data = example_data()
model = makemodel(data, :t, :id, :gub, :Y, [:X1, :X2], 
                  Dict(:X1 => false, :X2 => false), -15:-10, 1:5)
match!(model, wids=50)
balance!(model)
estimate!(model, dobayesfactor=false)

# Generate imputation results  
imputation = impute_results(model, model.matches, data, :t, :id)

# Comprehensive inspection
inspection = inspect_results(imputation, :Y, plot_percentage_changes=true)

# View plots
inspection.treatment_effects_plot    # Shows ATT over time
inspection.counterfactual_plot       # Shows causal story: observed vs counterfactual  
inspection.percentage_plot           # Shows effects as % changes

# Access summary statistics
summary = inspection.summary
println("Average treatment effect: ", round(summary.overall_att, digits=3))
println("Baseline balance: ", round(summary.baseline_difference, digits=3))

# Access underlying data for custom analysis  
results_df = inspection.results_df
custom_plot = plot_treatment_effects(results_df, title="My Custom ATT Plot")
```

# Notes
- This function coordinates the entire inspection workflow
- Data processing and plotting are separated for flexibility
- All plots contain proper labels and confidence intervals
- Use `plot_percentage_changes=true` for additional percentage change visualization
"""
function inspect_results(imputation_results::ImputationResults, outcome_var::Symbol;
                        plot_percentage_changes::Bool = false)
    
    # Validate inputs
    !isa(imputation_results, ImputationResults) && throw(ArgumentError("First argument must be ImputationResults"))
    outcome_var ∉ names(imputation_results.results) && throw(ArgumentError("Outcome variable $outcome_var not found in results"))
    
    # Data processing phase (separated from plotting) - now just sorts DataFrame
    results_df = prepare_treatment_effect_data(imputation_results)  # Same as others, just returns sorted DataFrame
    summary = calculate_overall_summary(imputation_results, outcome_var)
    
    # Plotting phase (pure visualization) - now takes DataFrames
    treatment_effects_plot = plot_treatment_effects(results_df)
    counterfactual_plot = plot_counterfactual_comparison(results_df, outcome_var)
    
    results = (
        treatment_effects_plot = treatment_effects_plot,
        counterfactual_plot = counterfactual_plot,
        summary = summary,
        results_df = results_df  # Give users access to the DataFrame
    )
    
    if plot_percentage_changes
        percentage_plot = plot_percentage_changes(results_df)  # Same DataFrame works for all!
        results = merge(results, (percentage_plot = percentage_plot,))
    end
    
    return results
end

"""
    create_inspection_dashboard(imputation_results::ImputationResults, outcome_var::Symbol)

Create a comprehensive dashboard with all inspection plots in a single figure.

# Arguments
- `imputation_results`: ImputationResults struct from impute_results()
- `outcome_var`: Symbol representing the outcome variable column name

# Returns
- `Figure`: Comprehensive dashboard with 4 panels:
  1. **Top Left**: Treatment effects (ATT) over time with confidence intervals
  2. **Top Right**: Observed vs counterfactual trajectories comparison  
  3. **Bottom Left**: Treatment effects as percentage changes
  4. **Bottom Right**: Summary statistics and balance diagnostics

# Example
```julia
# After running TSCSMethods analysis and generating imputation results
imputation = impute_results(model, model.matches, data, :t, :id)

# Create comprehensive dashboard
dashboard = create_inspection_dashboard(imputation, :Y)

# Save or display
save("analysis_dashboard.png", dashboard)
dashboard  # Display in notebook/REPL
```

# Dashboard Interpretation
- **ATT Panel**: Shows treatment effects over time. Look for consistent patterns and significant effects
- **Trajectory Panel**: Key causal story - distance between blue (observed) and black (counterfactual) lines shows treatment impact
- **Percentage Panel**: Treatment effects as % changes for easier interpretation of magnitude  
- **Summary Panel**: Key statistics including overall ATT and baseline balance assessment

# Notes
- Dashboard provides comprehensive view of results
- All panels use consistent time scales for easy comparison
- Confidence intervals help assess statistical uncertainty
- Summary statistics help evaluate matching quality and effect magnitude
"""
function create_inspection_dashboard(imputation_results::ImputationResults, outcome_var::Symbol)
    
    !PLOTTING_AVAILABLE && throw(ArgumentError("Dashboard functionality requires CairoMakie.jl"))
    
    # Validate inputs  
    !isa(imputation_results, ImputationResults) && throw(ArgumentError("First argument must be ImputationResults"))
    outcome_var ∉ names(imputation_results.results) && throw(ArgumentError("Outcome variable $outcome_var not found in results"))
    
    # Prepare data - now just one DataFrame for everything!
    results_df = prepare_treatment_effect_data(imputation_results)
    summary = calculate_overall_summary(imputation_results, outcome_var)
    
    # Create dashboard layout with enhanced styling
    fig = Figure(size = (1200, 800))
    
    # Top row: Treatment effects and counterfactual comparison
    ax1 = Axis(fig[1, 1], 
               title = "Treatment Effects Over Time", 
               xlabel = "Time Period", 
               ylabel = "ATT",
               xgridvisible = true,
               ygridvisible = true,
               xgridalpha = DEFAULT_PLOT_STYLE.grid_alpha,
               ygridalpha = DEFAULT_PLOT_STYLE.grid_alpha)
    
    ax2 = Axis(fig[1, 2], 
               title = "Observed vs Counterfactual",
               xlabel = "Time Period", 
               ylabel = string(outcome_var),
               xgridvisible = true,
               ygridvisible = true,
               xgridalpha = DEFAULT_PLOT_STYLE.grid_alpha,
               ygridalpha = DEFAULT_PLOT_STYLE.grid_alpha)
    
    # Bottom row: Percentage changes and summary
    ax3 = Axis(fig[2, 1], 
               title = "Treatment Effects (% Change)",
               xlabel = "Time Period", 
               ylabel = "Percent Change (%)",
               xgridvisible = true,
               ygridvisible = true,
               xgridalpha = DEFAULT_PLOT_STYLE.grid_alpha,
               ygridalpha = DEFAULT_PLOT_STYLE.grid_alpha)
    
    # Plot treatment effects with enhanced styling
    scatter!(ax1, results_df.f, results_df.att, 
             color = DEFAULT_PLOT_STYLE.colors.treatment, 
             markersize = DEFAULT_PLOT_STYLE.markersize)
    rangebars!(ax1, results_df.f, 
               results_df[!, Symbol("2.5%")], results_df[!, Symbol("97.5%")], 
               color = DEFAULT_PLOT_STYLE.colors.treatment, 
               linewidth = 2)
    hlines!(ax1, [0], 
            linestyle = :dash, 
            color = DEFAULT_PLOT_STYLE.colors.reference, 
            alpha = DEFAULT_PLOT_STYLE.reference_line_alpha,
            linewidth = 1.5)
    
    # Plot counterfactual comparison with enhanced styling
    # Confidence band first (behind lines)
    band!(ax2, results_df.f, 
          results_df.counterfactual_lower, results_df.counterfactual_upper,
          color = DEFAULT_PLOT_STYLE.colors.confidence, label = "95% CI")
    
    lines!(ax2, results_df.f, results_df[!, outcome_var], 
           color = DEFAULT_PLOT_STYLE.colors.observed, 
           linewidth = DEFAULT_PLOT_STYLE.linewidth, 
           label = "Observed")
    lines!(ax2, results_df.f, results_df.counterfactual_trajectory,
           color = DEFAULT_PLOT_STYLE.colors.counterfactual, 
           linewidth = DEFAULT_PLOT_STYLE.linewidth,
           linestyle = :dash,
           label = "Counterfactual")
    axislegend(ax2, position = :bottomright, framevisible = true, bgcolor = (:white, 0.9))
    
    # Plot percentage changes with significance highlighting
    significant = (results_df.pct_change_lower .> 0) .| (results_df.pct_change_upper .< 0)
    point_colors = [s ? :red : :orange for s in significant]
    
    scatter!(ax3, results_df.f, results_df.pct_change, 
             color = point_colors, 
             markersize = DEFAULT_PLOT_STYLE.markersize)
    rangebars!(ax3, results_df.f, 
               results_df.pct_change_lower, results_df.pct_change_upper, 
               color = point_colors,
               linewidth = 2)
    hlines!(ax3, [0], 
            linestyle = :dash, 
            color = DEFAULT_PLOT_STYLE.colors.reference, 
            alpha = DEFAULT_PLOT_STYLE.reference_line_alpha,
            linewidth = 1.5)
    
    # Add enhanced summary text with better formatting
    summary_text = """
    Summary Statistics
    
    Overall ATT: $(round(summary.overall_att, digits=3))
    
    Baseline Balance:
    • Difference: $(round(summary.baseline_difference, digits=3))
    • Treated Pre-avg: $(round(summary.treated_pretreatment_avg, digits=3))
    • Matched Pre-avg: $(round(summary.matched_pretreatment_avg, digits=3))
    
    Interpretation:
    • Baseline difference close to 0 indicates good matching
    • Overall ATT is the main causal estimate
    • Red points = statistically significant effects
    • Orange points = not statistically significant
    """
    
    Label(fig[2, 2], summary_text, 
          tellwidth = false, tellheight = false,
          halign = :left, valign = :top,
          fontsize = 11,
          lineheight = 1.3)
    
    return fig
end

## Generic Inspection Utilities

"""
    save_inspection_plots(inspection_results, output_dir::String = ".", prefix::String = "analysis")

Save all plots from inspect_results() to files with publication-quality settings.

# Arguments
- `inspection_results`: Output from inspect_results()
- `output_dir`: Directory to save plots (default: current directory)
- `prefix`: File prefix for saved plots (default: "analysis")

# Example
```julia
# Run inspection
inspection = inspect_results(imputation, :Y, plot_percentage_changes=true)

# Save all plots
save_inspection_plots(inspection, "results/", "covid_analysis")
# Creates: covid_analysis_treatment_effects.png, covid_analysis_counterfactual.png, etc.
```
"""
function save_inspection_plots(inspection_results, output_dir::String = ".", prefix::String = "analysis")
    !PLOTTING_AVAILABLE && throw(ArgumentError("Saving plots requires CairoMakie.jl"))
    
    # Ensure output directory exists
    !isdir(output_dir) && mkpath(output_dir)
    
    # Save treatment effects plot
    if haskey(inspection_results, :treatment_effects_plot)
        save(joinpath(output_dir, "$(prefix)_treatment_effects.png"), 
             inspection_results.treatment_effects_plot, px_per_unit = 2)
        println("Saved: $(prefix)_treatment_effects.png")
    end
    
    # Save counterfactual plot
    if haskey(inspection_results, :counterfactual_plot)
        save(joinpath(output_dir, "$(prefix)_counterfactual.png"), 
             inspection_results.counterfactual_plot, px_per_unit = 2)
        println("Saved: $(prefix)_counterfactual.png")
    end
    
    # Save percentage plot if it exists
    if haskey(inspection_results, :percentage_plot)
        save(joinpath(output_dir, "$(prefix)_percentage_changes.png"), 
             inspection_results.percentage_plot, px_per_unit = 2)
        println("Saved: $(prefix)_percentage_changes.png")
    end
    
    return nothing
end

"""
    export_results_csv(inspection_results, output_path::String)

Export the results DataFrame to CSV for further analysis in other software.

# Arguments
- `inspection_results`: Output from inspect_results()
- `output_path`: Path for the CSV file

# Example
```julia
inspection = inspect_results(imputation, :Y)
export_results_csv(inspection, "analysis_results.csv")
```
"""
function export_results_csv(inspection_results, output_path::String)
    if haskey(inspection_results, :results_df)
        CSV.write(output_path, inspection_results.results_df)
        println("Results exported to: $output_path")
    else
        throw(ArgumentError("No results_df found in inspection_results"))
    end
    return nothing
end

"""
    quick_inspection(imputation_results::ImputationResults, outcome_var::Symbol)

Quick inspection function that prints summary statistics and optionally shows plots.

# Arguments
- `imputation_results`: Output from impute_results()
- `outcome_var`: Symbol for outcome variable

# Returns
- Summary statistics printed to console
- Basic diagnostic information about matching quality

# Example
```julia
imputation = impute_results(model, matches, data, :t, :id)
quick_inspection(imputation, :Y)
```
"""
function quick_inspection(imputation_results::ImputationResults, outcome_var::Symbol)
    summary = calculate_overall_summary(imputation_results, outcome_var)
    results_df = prepare_treatment_effect_data(imputation_results)
    
    println("=== TSCSMethods Quick Inspection ===")
    println()
    
    println("Sample Information:")
    println("  Time periods: $(minimum(results_df.f)) to $(maximum(results_df.f))")
    println("  Number of periods: $(nrow(results_df))")
    println()
    
    println("Treatment Effects:")
    println("  Overall ATT: $(round(summary.overall_att, digits=4))")
    println("  ATT Range: $(round(minimum(results_df.att), digits=3)) to $(round(maximum(results_df.att), digits=3))")
    
    # Check for significant effects
    significant_periods = sum((results_df[!, Symbol("2.5%")] .> 0) .| (results_df[!, Symbol("97.5%")] .< 0))
    println("  Significant periods: $significant_periods / $(nrow(results_df))")
    println()
    
    println("Matching Quality:")
    println("  Baseline difference: $(round(summary.baseline_difference, digits=4))")
    println("  Treated pre-treatment avg: $(round(summary.treated_pretreatment_avg, digits=3))")
    println("  Matched pre-treatment avg: $(round(summary.matched_pretreatment_avg, digits=3))")
    
    # Assess matching quality
    baseline_quality = abs(summary.baseline_difference) < 0.1 ? "Good" : 
                      abs(summary.baseline_difference) < 0.5 ? "Moderate" : "Poor"
    println("  Baseline balance quality: $baseline_quality")
    println()
    
    if PLOTTING_AVAILABLE
        println("For detailed plots, use:")
        println("  inspection = inspect_results(imputation, :$outcome_var)")
        println("  inspection.treatment_effects_plot")
        println("  inspection.counterfactual_plot")
    else
        println("Install CairoMakie.jl for plotting capabilities")
    end
    
    return summary
end

"""
    compare_strata(imputation_results_list::Vector, outcome_var::Symbol, stratum_names::Vector{String})

Compare results across multiple strata or model specifications.

# Arguments
- `imputation_results_list`: Vector of ImputationResults from different strata
- `outcome_var`: Symbol for outcome variable
- `stratum_names`: Names for each stratum for labeling

# Returns
- Comparison plot showing ATT estimates across strata
- Summary table of key statistics

# Example
```julia
# Run analysis for each stratum
imputation1 = impute_results(model, matches, data, :t, :id, stratum=1)
imputation2 = impute_results(model, matches, data, :t, :id, stratum=2)

# Compare results
comparison = compare_strata([imputation1, imputation2], :Y, ["Urban", "Rural"])
```
"""
function compare_strata(imputation_results_list::Vector, outcome_var::Symbol, stratum_names::Vector{String})
    !PLOTTING_AVAILABLE && throw(ArgumentError("Stratum comparison plotting requires CairoMakie.jl"))
    
    length(imputation_results_list) != length(stratum_names) && 
        throw(ArgumentError("Number of imputation results must match number of stratum names"))
    
    # Prepare data for each stratum
    all_summaries = []
    all_results = []
    
    for (i, imputation) in enumerate(imputation_results_list)
        summary = calculate_overall_summary(imputation, outcome_var)
        results_df = prepare_treatment_effect_data(imputation)
        results_df[!, :stratum] .= stratum_names[i]
        
        push!(all_summaries, summary)
        push!(all_results, results_df)
    end
    
    # Combine all results
    combined_df = vcat(all_results...)
    
    # Create comparison plot
    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1, 1], 
              title = "Treatment Effects Comparison Across Strata",
              xlabel = "Time Period", 
              ylabel = "ATT",
              xgridvisible = true,
              ygridvisible = true)
    
    colors = [:blue, :red, :green, :orange, :purple]
    
    for (i, stratum_name) in enumerate(stratum_names)
        stratum_data = filter(row -> row.stratum == stratum_name, combined_df)
        color = colors[mod1(i, length(colors))]
        
        scatter!(ax, stratum_data.f, stratum_data.att,
                 color = color, markersize = 6, label = stratum_name)
        rangebars!(ax, stratum_data.f, 
                   stratum_data[!, Symbol("2.5%")], stratum_data[!, Symbol("97.5%")],
                   color = color, alpha = 0.7)
    end
    
    hlines!(ax, [0], linestyle = :dash, color = :gray, alpha = 0.5)
    axislegend(ax)
    
    # Print summary comparison
    println("=== Stratum Comparison Summary ===")
    for (i, (name, summary)) in enumerate(zip(stratum_names, all_summaries))
        println("$name:")
        println("  Overall ATT: $(round(summary.overall_att, digits=4))")
        println("  Baseline difference: $(round(summary.baseline_difference, digits=4))")
    end
    
    return (plot = fig, summaries = all_summaries, data = combined_df)
end

## Utility Functions

"""
    change_pct(val, attval)

Calculate percentage change from baseline value.

# Arguments
- `val`: Baseline value
- `attval`: Treatment effect (difference from baseline)

# Returns  
- `Float64`: Percentage change (100 * attval / val)

# Example
```julia
baseline = 100.0
treatment_effect = 10.0
pct_change = change_pct(baseline, treatment_effect)  # Returns 10.0 (10% increase)
```
"""
function change_pct(val, attval)
    return 100 * attval / val
end

#=
## Complete Workflow Example

Here's a comprehensive example showing the full TSCSMethods + inspection workflow:

```julia
using TSCSMethods

# 1. Load and prepare data
data = example_data()

# 2. Create and fit TSCSMethods model  
model = makemodel(
    data,                    # Panel data
    :t,                      # Time variable  
    :id,                     # Unit ID variable
    :gub,                    # Treatment variable
    :Y,                      # Outcome variable
    [:X1, :X2],             # Covariates
    Dict(:X1 => false, :X2 => false),  # Time-varying specification
    -15:-10,                 # Pre-treatment window
    1:5                      # Post-treatment window
)

# 3. Run matching, balancing, and estimation
match!(model, wids=50)              # Match treated units to controls
balance!(model)                      # Balance covariates
estimate!(model, dobayesfactor=false)  # Estimate treatment effects

# 4. Generate counterfactual analysis
imputation = impute_results(model, model.matches, data, :t, :id)

# 5. Comprehensive inspection (Method 1: Individual plots)
inspection = inspect_results(imputation, :Y, plot_percentage_changes=true)

# View individual plots
inspection.treatment_effects_plot  # ATT over time
inspection.counterfactual_plot     # Observed vs counterfactual  
inspection.percentage_plot         # Effects as percentages

# Access summary statistics
println("Overall ATT: ", round(inspection.summary.overall_att, digits=3))
println("Baseline balance: ", round(inspection.summary.baseline_difference, digits=3))

# 6. Comprehensive dashboard (Method 2: All-in-one view)  
dashboard = create_inspection_dashboard(imputation, :Y)
save("analysis_results.png", dashboard)

# 7. Custom analysis using prepared data
results_df = prepare_treatment_effect_data(imputation)  # Just returns sorted DataFrame

# Create custom visualizations
custom_plot = plot_treatment_effects(
    results_df,
    title="My Custom Treatment Effects",
    ylabel="Change in Outcome",
    color=:red
)
```

## Interpretation Guide

### Treatment Effects Plot
- **Positive values**: Treatment increased the outcome
- **Negative values**: Treatment decreased the outcome  
- **Confidence intervals**: Wider intervals indicate more uncertainty
- **Time patterns**: Look for immediate vs. delayed effects

### Counterfactual Comparison Plot
- **Blue line (Observed)**: What actually happened to treated units
- **Black line (Counterfactual)**: What would have happened without treatment
- **Gray band**: Uncertainty around counterfactual estimates
- **Gap between lines**: The causal effect of treatment

### Summary Statistics
- **Overall ATT**: Average treatment effect across all time periods
- **Baseline difference**: Pre-treatment difference between treated and matched units
  - Should be close to 0 for good matching
- **Pre-treatment averages**: Help assess matching quality

### Percentage Changes
- Easier interpretation of treatment magnitude
- Shows treatment effects as % of baseline outcome levels
- Useful for communicating results to non-technical audiences

## Best Practices

1. **Always check baseline balance** (`baseline_difference`) before interpreting results
2. **Look for pre-treatment trends** in the counterfactual comparison plot
3. **Assess statistical significance** using confidence intervals
4. **Consider practical significance** alongside statistical significance
5. **Use percentage changes** for easier interpretation and communication
6. **Save high-resolution plots** using `save("filename.png", plot, px_per_unit=2)`

## Troubleshooting

- **"Plotting functionality requires CairoMakie.jl"**: Install with `Pkg.add("CairoMakie")`
- **"Outcome variable not found"**: Check that the Symbol matches a column name exactly
- **Missing confidence intervals**: Ensure `estimate!()` was run with proper bootstrap settings
- **Poor baseline balance**: Consider different matching specifications or caliper restrictions
- **Jagged counterfactual**: May indicate insufficient matched controls; increase `wids` in `match!()`

=#

# Note: The inspection.jl file has been completely modernized with new functions:
# - Use inspect_results() for comprehensive analysis
# - Use quick_inspection() for command-line summaries  
# - Use create_inspection_dashboard() for publication-ready dashboards
# - See function documentation above for detailed examples

## Pure Plotting Functions (work with prepared data)

"""
    plot_treatment_effects(results_df; kwargs...)

Plot treatment effects over time with confidence intervals.
Takes DataFrame from prepare_treatment_effect_data().

# Arguments
- `results_df`: DataFrame with columns `f`, `att`, `2.5%`, `97.5%`
- `title`: Plot title (default: "Treatment Effects Over Time")
- `xlabel`: X-axis label (default: "Time Period") 
- `ylabel`: Y-axis label (default: "ATT")
- `markersize`: Size of scatter points (default: from DEFAULT_PLOT_STYLE)
- `color`: Color for points and error bars (default: steelblue)
- `show_grid`: Whether to show grid lines (default: true)
- `show_zero_line`: Whether to show reference line at zero (default: true)

# Returns
- `Figure`: Makie figure object
"""

# Plotting defaults
const DEFAULT_PLOT_STYLE = (
    markersize = 8,
    linewidth = 2.5,
    ci_alpha = 0.3,
    grid_alpha = 0.2,
    reference_line_alpha = 0.6,
    colors = (
        treatment = :steelblue,
        counterfactual = :black,
        observed = :darkblue,
        confidence = (:gray, 0.3),
        reference = :gray
    )
)

function plot_treatment_effects(results_df::DataFrame; 
                               title::String = "Treatment Effects Over Time",
                               xlabel::String = "Time Period",
                               ylabel::String = "ATT",
                               markersize::Int = DEFAULT_PLOT_STYLE.markersize,
                               color = DEFAULT_PLOT_STYLE.colors.treatment,
                               show_grid::Bool = true,
                               show_zero_line::Bool = true)
    
    !PLOTTING_AVAILABLE && throw(ArgumentError("Plotting functionality requires CairoMakie.jl"))
    
    # Input validation - check required columns
    required_cols = [:f, :att, Symbol("2.5%"), Symbol("97.5%")]
    for col in required_cols
        col ∉ names(results_df) && throw(ArgumentError("DataFrame must have column: $col"))
    end
    
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], 
              title = title, 
              xlabel = xlabel, 
              ylabel = ylabel,
              xgridvisible = show_grid,
              ygridvisible = show_grid,
              xgridalpha = DEFAULT_PLOT_STYLE.grid_alpha,
              ygridalpha = DEFAULT_PLOT_STYLE.grid_alpha)
    
    # Plot ATT estimates with enhanced styling
    scatter!(ax, results_df.f, results_df.att, 
             color = color, 
             markersize = markersize,
             label = "Point Estimate")
    
    # Enhanced confidence intervals
    rangebars!(ax, results_df.f, 
               results_df[!, Symbol("2.5%")], results_df[!, Symbol("97.5%")],
               color = color,
               linewidth = 2,
               label = "95% Confidence Interval")
    
    # Optional reference line at zero
    if show_zero_line
        hlines!(ax, [0], 
                linestyle = :dash, 
                color = DEFAULT_PLOT_STYLE.colors.reference, 
                alpha = DEFAULT_PLOT_STYLE.reference_line_alpha,
                linewidth = 1.5)
    end
    
    # Add legend
    axislegend(ax, position = :topright)
    
    return fig
end

"""
    plot_counterfactual_comparison(results_df, outcome_var; kwargs...)

Plot observed vs counterfactual trajectories.
Takes DataFrame from prepare_counterfactual_comparison_data().

# Arguments
- `results_df`: DataFrame with columns `f`, `[outcome_var]`, `counterfactual_trajectory`, etc.
- `outcome_var`: Symbol for outcome variable column name
- `title`: Plot title (default: "Observed vs Counterfactual Trajectories")
- `xlabel`: X-axis label (default: "Time Period")
- `ylabel`: Y-axis label (default: string of outcome_var)
"""
function plot_counterfactual_comparison(results_df::DataFrame, outcome_var::Symbol;
                                      title::String = "Observed vs Counterfactual Trajectories",
                                      xlabel::String = "Time Period",
                                      ylabel::String = string(outcome_var),
                                      show_grid::Bool = true,
                                      add_vertical_treatment_line::Bool = false,
                                      treatment_start::Union{Nothing, Real} = nothing)
    
    !PLOTTING_AVAILABLE && throw(ArgumentError("Plotting functionality requires CairoMakie.jl"))
    
    # Input validation - check required columns
    required_cols = [:f, outcome_var, :counterfactual_trajectory, :counterfactual_lower, :counterfactual_upper]
    for col in required_cols
        col ∉ names(results_df) && throw(ArgumentError("DataFrame must have column: $col"))
    end
    
    fig = Figure(size = (900, 600))
    ax = Axis(fig[1, 1], 
              title = title, 
              xlabel = xlabel, 
              ylabel = ylabel,
              xgridvisible = show_grid,
              ygridvisible = show_grid,
              xgridalpha = DEFAULT_PLOT_STYLE.grid_alpha,
              ygridalpha = DEFAULT_PLOT_STYLE.grid_alpha)
    
    # Enhanced confidence band first (so it's behind lines)
    band!(ax, results_df.f, 
          results_df.counterfactual_lower,
          results_df.counterfactual_upper,
          color = DEFAULT_PLOT_STYLE.colors.confidence, 
          label = "95% CI")
    
    # Plot observed trajectory with enhanced styling
    lines!(ax, results_df.f, results_df[!, outcome_var], 
           color = DEFAULT_PLOT_STYLE.colors.observed, 
           linewidth = DEFAULT_PLOT_STYLE.linewidth, 
           label = "Observed (Treated)")
           
    # Plot counterfactual trajectory with enhanced styling
    lines!(ax, results_df.f, results_df.counterfactual_trajectory,
           color = DEFAULT_PLOT_STYLE.colors.counterfactual, 
           linewidth = DEFAULT_PLOT_STYLE.linewidth,
           linestyle = :dash,
           label = "Counterfactual")
    
    # Optional vertical line at treatment start
    if add_vertical_treatment_line && !isnothing(treatment_start)
        vlines!(ax, [treatment_start], 
                linestyle = :dot, 
                color = :red, 
                alpha = 0.7,
                linewidth = 2,
                label = "Treatment Start")
    end
    
    # Enhanced legend with better positioning
    axislegend(ax, position = :bottomright, framevisible = true, 
               bgcolor = (:white, 0.9))
    
    return fig
end

"""
    plot_percentage_changes(results_df; kwargs...)

Plot treatment effects as percentage changes.
Takes DataFrame from prepare_percentage_change_data().

# Arguments  
- `results_df`: DataFrame with columns `f`, `pct_change`, `pct_change_lower`, `pct_change_upper`
- `title`: Plot title (default: "Treatment Effects (% Change)")
- `xlabel`: X-axis label (default: "Time Period")
- `ylabel`: Y-axis label (default: "Percent Change")  
- `markersize`: Size of scatter points (default: 8)
"""
function plot_percentage_changes(results_df::DataFrame;
                                title::String = "Treatment Effects as Percentage Changes",
                                xlabel::String = "Time Period", 
                                ylabel::String = "Percent Change (%)",
                                markersize::Int = DEFAULT_PLOT_STYLE.markersize,
                                show_grid::Bool = true,
                                show_zero_line::Bool = true,
                                highlight_significant::Bool = true)
    
    !PLOTTING_AVAILABLE && throw(ArgumentError("Plotting functionality requires CairoMakie.jl"))
    
    # Input validation - check required columns
    required_cols = [:f, :pct_change, :pct_change_lower, :pct_change_upper]
    for col in required_cols
        col ∉ names(results_df) && throw(ArgumentError("DataFrame must have column: $col"))
    end
    
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], 
              title = title, 
              xlabel = xlabel, 
              ylabel = ylabel,
              xgridvisible = show_grid,
              ygridvisible = show_grid,
              xgridalpha = DEFAULT_PLOT_STYLE.grid_alpha,
              ygridalpha = DEFAULT_PLOT_STYLE.grid_alpha)
    
    # Determine point colors based on significance if requested
    point_colors = if highlight_significant
        # Check if confidence intervals exclude zero
        significant = (results_df.pct_change_lower .> 0) .| (results_df.pct_change_upper .< 0)
        [s ? :red : :orange for s in significant]
    else
        fill(:red, nrow(results_df))
    end
    
    # Plot percentage changes with enhanced styling
    scatter!(ax, results_df.f, results_df.pct_change,
             color = point_colors, 
             markersize = markersize,
             label = highlight_significant ? "● Significant  ● Not Significant" : "Point Estimate")
    
    # Enhanced confidence intervals
    rangebars!(ax, results_df.f, 
               results_df.pct_change_lower, results_df.pct_change_upper,
               color = point_colors,
               linewidth = 2)
    
    # Optional reference line at zero
    if show_zero_line
        hlines!(ax, [0], 
                linestyle = :dash, 
                color = DEFAULT_PLOT_STYLE.colors.reference, 
                alpha = DEFAULT_PLOT_STYLE.reference_line_alpha,
                linewidth = 1.5,
                label = "No Effect")
    end
    
    # Add legend only if showing significance colors or zero line
    if highlight_significant || show_zero_line
        axislegend(ax, position = :topright)
    end
    
    return fig
end
