# inspection.jl
#
# CURRENT STATUS: This file contains legacy code that needs major refactoring
# 
# DEPRECATED FUNCTIONS REMOVED:
# - pretreatment(): Used hardcoded 'fips'/'running' columns and magic numbers
# - inspection(): Referenced undefined variables and used old API
# - plot_inspection(): Used incompatible column names
#
# INCOMPATIBLE FUNCTIONS:
# - figure_6(): Expects old column names from before imputation.jl refactor
# - fill_att!(): Uses old confidence interval column names  
#
# WORKING UTILITY FUNCTIONS:
# - change_pct(): Simple percentage calculation (works)
# - fill_axis!(), fill_baxis!(): Basic plotting helpers (work with correct data)
# - model_balance_plot(): Balance plotting (works but has hardcoded styling)

function change_pct(val, attval)
    return 100 * attval / val
end

# DEPRECATED: Broken inspection functions with undefined variables
# TODO: Replace with properly designed inspection functions
# 
# Issues with old functions:
# - References undefined 'dat', 'objet' variables
# - Calls non-existent 'load_object' function  
# - Uses old impute_results signature
# - Returns inconsistent values

# DEPRECATED: Dataset-specific function with hardcoded column names
# TODO: Replace with generic pretreatment analysis function
# function pretreatment(matches, oc, dat, time_col, unit_col)
#     This function used hardcoded 'fips' and 'running' column names
#     and magic number 20 for time periods. Needs complete rewrite.

## plotting

function fill_axis!(ax1, t, dm, variables)
    lines!(ax1, t, dm[!, variables[1]], color = :blue)
    lines!(ax1, t, dm[!, variables[2]], color = :black)
    band!(t, dm[!, variables[3]], dm[!, variables[4]], color = (:grey, 0.3)) 
end

function fill_att!(ax2, t, dm; msize = 8)
    scatter!(ax2, t, dm[!, :att], color = :black, markersize= msize)
    rangebars!(ax2, t, dm[!, Symbol("2.5%")], dm[!, Symbol("97.5%")])
end

# DEPRECATED: Old plotting function with hardcoded column names  
# TODO: Replace with functions that work with ImputationResults struct
# function plot_inspection(imputation_results::ImputationResults, ...)

# DEPRECATED: Dataset-specific plotting function
# TODO: Modernize to work with new ImputationResults and remove hardcoded assumptions
# INCOMPATIBLE: This function expects old column names that no longer exist:
#   :counter_post -> :counterfactual_trajectory
#   :counter_post_lwr -> :counterfactual_lower  
#   :counter_post_upr -> :counterfactual_upper
#   :pct -> :pct_change
#   :pct_lo -> :pct_change_lower
#   :pct_hi -> :pct_change_upper
function figure_6(dm, overallestimate, oc, stratum; msize = 8, plot_pct = false)
    fig = Figure();
    oe = overallestimate[stratum];

    treated_observed_mean = mean(dm[!, oc]);
    overall_att = oe[1];
    lo_att = oe[2][1]; hi_att = oe[2][3];

    # BROKEN: These column names no longer exist after imputation.jl refactoring
    variables = [oc, :counter_post, :counter_post_lwr, :counter_post_upr];
    t = sort(unique(dm.f));

    ftop = fig[1, 1:2] = GridLayout();
    fbot = fig[2, 1:2] = GridLayout();
    if plot_pct
        fbotbot = fig[3,1:2] = GridLayout();
    end

    ax1 = Axis(fbot[1,1], ylabel = string(oc), xlabel = "Day");
    fill_axis!(ax1, t, dm, variables)
    mn, mx = extrema(t)
    xlims!(ax1, [mn - 0.75, mx + 0.25])

    ax2 = Axis(ftop[1,1], ylabel = string(oc));
    fill_att!(ax2, t, dm)
    xlims!(ax2, [mn - 0.75, mx + 0.25])

    ax3 = Axis(ftop[1,2]);

    scatter!(ax3, [9], [overall_att], color = :black, markersize = msize)
    rangebars!(ax3, [9], lo_att, hi_att)

    linkyaxes!(ax2, ax3)
    
    hideydecorations!(ax3, grid = false)
    hidexdecorations!(ax3)
    
    ## overall

    # mean
    # mean(treated outcome values) - mean(att estimate)
    # att > 0 -> treated values would be higher
    # att < 0 -> treated values would be lower
    br = treated_observed_mean - overall_att
    
    # 95% CI
    # mean(treated outcome value) - [lwr att est, upr att est]
    br_lo = treated_observed_mean - lo_att
    br_hi = treated_observed_mean - hi_att

    ax4 = Axis(fbot[1,2]);

    # blue: average for the treated counties
    # black: counterfactual average if they had not been tretaed
    scatter!(
        ax4, [3.25, 2.75], [br, treated_observed_mean],
        color = [:black, :blue], markersize = msize
    )
    # add counterfactual at CI bounds
    rangebars!(ax4, [3.25], br_lo, br_hi)
    xlims!(ax4, [2.5, 3.5])

    hideydecorations!(ax4, grid = false)
    hidexdecorations!(ax4)

    linkyaxes!(ax1, ax4)

    colsize!(ftop, 2, Auto(0.04))
    colsize!(fbot, 2, Auto(0.04))

    for (label, layout) in zip(["a", "b"], [ftop, fbot])
        Label(layout[1, 1, TopLeft()], label,
            textsize = 26,
            padding = (0, 5, 5, 0),
            halign = :right)
    end


    if plot_pct
        ax5 = Axis(fbotbot[1,1], ylabel = "Pct. difference");
        # BROKEN: :pct should be :pct_change
        scatter!(ax5, t, dm[!, :pct], color = :black, markersize = msize)
        # BROKEN: :pct_lo/:pct_hi should be :pct_change_lower/:pct_change_upper
        rangebars!(ax5, t, dm[!, :pct_lo], dm[!, :pct_hi])

        ax6 = Axis(fbotbot[1,2]);
        # BROKEN: dm.pct should be dm.pct_change  
        cm = mean(dm.pct) # == change_pct(treated_observed_mean, oe[1])
        clo = change_pct(treated_observed_mean, oe[2][1])
        chi = change_pct(treated_observed_mean, oe[2][2])

        scatter!(
            ax6, [9], [cm],
            color = [:black], markersize = msize
        )
        # add counterfactual at CI bounds
        rangebars!(ax6, [9], clo, chi)
        # xlims!(ax6, [2.5, 3.5])

        linkyaxes!(ax5, ax6)
    
        hideydecorations!(ax6, grid = false)
        hidexdecorations!(ax6)

        colsize!(fbotbot, 2, Auto(0.04))
    end

    return fig
end

function fill_baxis!(bax, d_gb, t)
    for (k, v) in sort(d_gb[1])
        bl = if typeof(v) <: Vector
            v
        else
            fill(v, length(t))
        end
        lines!(bax, bl, label = string(k))
    end
end

function model_balance_plot(d_gb)

    bfig = Figure();
    bax = Axis(bfig[1,1], xlabel = "Day", ylabel = "Balance score");

    fill_baxis!(bax, d_gb, 10:40)
    Legend(bfig[1,2], bax, "Covariates", framevisible = false);
    ylims!(bax, -0.15, 0.15)

    hlines!(bax, [-0.1, 0.1]; linestyle = :dash, color = :black)
    hideydecorations!(bax, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(bax, ticks = false, ticklabels = false, label = false)

    bax.xticks = (1:5:30, string.(collect(-30:5:-1)))

    return bfig
end
