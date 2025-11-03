##### plots #####
# requires results from the main file

using CairoMakie
include("functions_plots.jl")

##### comparison plot for Roulstone vs. Rothwell scaling #####

# α range
x = 0:0.01:1;
# β Roulstone
y = x;
# β Rothwell
z = 1 .+ log.(x)/log(2);

# define plot
fig_theory = Figure();

ax_theory = Axis(fig_theory[1,1], xlabel = "α", ylabel = "β(α)");
xlims!(0, 1)
ylims!(-2, 1)

roulstone = lines!(x, y, label = "Roulstone", linewidth = 3, color = :darkblue);
rothwell = lines!(x, z, label = "Rothwell", linewidth = 3, color = :green);
hlines!(ax_theory, [0], color = :gray);

Legend(fig_theory[1, 1],
    [roulstone, rothwell],
    [L"β^\text{Roulstone}", L"β^\text{Rothwell}"],
    tellheight = false,
    tellwidth = false,
    halign = :right, valign = :bottom,
    framevisible = false, orientation = :vertical);

fig_theory
save("$outputpath/fig-theory.pdf", fig_theory);

##### comparison plot for investment cost from manufacturers vs. estimation #####

# choose scaling parameters for the plot
scaling_plot = [0.20, 0.75];

# Original plot (all reactors)
# fig_invest_comparison = investment_plot(pjs, scaling_plot)
# save("$outputpath/fig-investment_comparison.pdf", fig_invest_comparison);

# New: Grouped by scale (Micro/SMR/Large)
# Use Base.invokelatest to avoid Julia 1.12 world age issues
fig_invest_comparison_by_scale = Base.invokelatest(investment_plot_by_scale, pjs, scaling_plot)
save("$outputpath/fig-investment_comparison_by_scale.pdf", fig_invest_comparison_by_scale);

##### histogram plots for comparison of estimation approaches #####

hist_invest = Figure();

for i in 1:3, j in 1:5
    hist_invest_plot(n, wacc, electricity_price_mean, pjs[j+5*(i-1)], i, j, hist_invest)
end

Legend(hist_invest[4,1:5],
    [roulstone, rothwell],
    ["Roulstone", "Rothwell"],
    framevisible = false, orientation = :horizontal)

hist_invest
save("$outputpath/fig-histogram_investment.pdf", hist_invest);

##### probability density plots for comparison of estimation approaches #####

density_invest = Figure();

for i in 1:3, j in 1:5
    density_invest_plot(n, wacc, electricity_price_mean, pjs[j+5*(i-1)], i, j, density_invest)
end

Legend(density_invest[4,1:5],
    [roulstone, rothwell],
    ["Roulstone", "Rothwell"],
    framevisible = false, orientation = :horizontal)

density_invest
save("$outputpath/fig-density_investment.pdf", density_invest);

##### boxplots Monte Carlo simulation results #####
# requires results for all 15 reactor concepts

fig_mcs_npv = mcs_plot(npv_results, "NPV", "[USD/MW]")
fig_mcs_lcoe = mcs_plot(lcoe_results, "LCOE", "[USD/MWh]")

save("$outputpath/fig-mcs_npv-$opt_scaling.pdf", fig_mcs_npv);
save("$outputpath/fig-mcs_lcoe-$opt_scaling.pdf", fig_mcs_lcoe);

##### heatmaps sensitivity indices #####
# requires sensitvity results

# Original plots (all reactors - may be too crowded)
# fig_si_npv = si_plot(si_npv_results, "NPV")
# fig_si_lcoe = si_plot(si_lcoe_results, "LCOE")
# save("$outputpath/fig-si_npv-$opt_scaling.pdf", fig_si_npv);
# save("$outputpath/fig-si_lcoe-$opt_scaling.pdf", fig_si_lcoe);

# New: Grouped by scale (Micro/SMR/Large)
# Use Base.invokelatest to avoid Julia 1.12 world age issues
fig_si_npv_by_scale = Base.invokelatest(si_plot_by_scale, si_npv_results, "NPV Sensitivity Indices", pjs)
fig_si_lcoe_by_scale = Base.invokelatest(si_plot_by_scale, si_lcoe_results, "LCOE Sensitivity Indices", pjs)

save("$outputpath/fig-si_npv_by_scale-$opt_scaling.pdf", fig_si_npv_by_scale);
save("$outputpath/fig-si_lcoe_by_scale-$opt_scaling.pdf", fig_si_lcoe_by_scale);

##### lcoe comparison plot #####
# requires results for all 15 reactor concepts

using CSV, DataFrames

# read LCOE data from CSV into a dataframe
lcoe_dat = CSV.File("$inputpath/lcoe_data.csv") |> DataFrame;

# read LCOE data from project simulations
lcoe_bounds = select(lcoe_summary, [:q25, :q75]);

# collect plot data
lcoe_plot_data = vcat(
    select(lcoe_dat, [:technology, :lower_bound, :upper_bound]),
    DataFrame(technology = ["BWR & PWR SMRs", "HTR SMRs", "SFR SMRs"],
    lower_bound = [minimum(lcoe_bounds[1:9,:q25]), minimum(lcoe_bounds[10:12,:q25]), minimum(lcoe_bounds[13:15,:q25])],
    upper_bound = [maximum(lcoe_bounds[1:9,:q75]), maximum(lcoe_bounds[10:12,:q75]), maximum(lcoe_bounds[13:15,:q75])])
);

# define LCOE plot
if opt_scaling == "manufacturer"
    plot_scaling = "Manufacturer";
elseif opt_scaling == "roulstone"
    plot_scaling = "Roulstone";
elseif opt_scaling == "rothwell"
    plot_scaling = "Rothwell";
elseif opt_scaling == "uniform"
    plot_scaling = "uniform";
else
    @error("scaling not defined")
end
xlabel = "[USD/MWh]";
yticks = lcoe_plot_data[!,:technology];

col = vcat(fill(1,8),fill(2,4),3,4,5);
colormap = [:darkgreen, :darkblue];

fig_lcoe_comparison = Figure();
ax_lcoe = Axis(fig_lcoe_comparison[1,1], yticks = (1:length(yticks), yticks), xscale = log10, xlabel = xlabel);

xlims!(10, 25000)

rangebars!(ax_lcoe, 1:length(yticks), lcoe_plot_data[!,2], lcoe_plot_data[!,3], linewidth = 6, whiskerwidth = 8, direction = :x, color = col);
hlines!(ax_lcoe, [8.5, 12.5], linestyle = :dash, color = :red);
text!([15000,15000,16], [4, 10.5, 14]; text = ["Renewables\n(LAZARD)", "Conventionals\n(LAZARD)", "SMR Tech.\n($plot_scaling)"], align = (:center, :center), justification = :center, rotation = π/2);

text!(lcoe_plot_data[!,2], 1:length(yticks), text = string.(round.(Int,lcoe_plot_data[!,2])), align = (:right, :center), offset = (-10,0));
text!(lcoe_plot_data[!,3], 1:length(yticks), text = string.(round.(Int,lcoe_plot_data[!,3])), align = (:left, :center), offset = (10,0));

Label(fig_lcoe_comparison[1, 1, Top()], "LCOE Comparison", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

fig_lcoe_comparison
save("$outputpath/fig-lcoe_comparison-$opt_scaling.pdf", fig_lcoe_comparison);
##### LCOE histogram grouped by scale (Micro/SMR/Large) #####
# requires lcoe_results from simulation and pjs vector

fig_lcoe_scale_hist = lcoe_scale_histogram(lcoe_results, pjs)
save("$outputpath/fig-lcoe_scale_histogram-$opt_scaling.pdf", fig_lcoe_scale_hist);

##### LCOE threshold probability plot #####
# Shows cumulative probability: P(LCOE ≤ threshold) by scale
# Use Base.invokelatest to avoid Julia 1.12 world age issues

# Define thresholds (0 to 300 USD/MWh in steps of 20)
lcoe_thresholds = collect(0.0:20.0:300.0)  # Float64 values required
fig_lcoe_threshold_prob = Base.invokelatest(lcoe_threshold_probability_plot, lcoe_results, pjs; thresholds=lcoe_thresholds)
save("$outputpath/fig-lcoe_threshold_probability-$opt_scaling.pdf", fig_lcoe_threshold_prob);

##### Regional comparison plots #####
# Compare LCOE distributions across regions
# Shows both full dataset and Large reactor specific comparisons

@info("Generating regional comparison plots")

# ALL SCALES: Show complete dataset grouped by region
@info("Generating regional plots for all scales")
fig_regional_all = Base.invokelatest(mcs_plot_regional, lcoe_results, pjs; scale_filter="All")
save("$outputpath/fig-regional_comparison_all-$opt_scaling.pdf", fig_regional_all);

fig_regional_combined_all = Base.invokelatest(mcs_plot_regional_combined, lcoe_results, pjs; scale_filter="All")
save("$outputpath/fig-regional_combined_all-$opt_scaling.pdf", fig_regional_combined_all);

# LARGE ONLY: Specific comparison for deployment learning analysis
@info("Generating regional plots for Large reactors only")
fig_regional_large = Base.invokelatest(mcs_plot_regional, lcoe_results, pjs; scale_filter="Large")
save("$outputpath/fig-regional_comparison_large-$opt_scaling.pdf", fig_regional_large);

fig_regional_combined_large = Base.invokelatest(mcs_plot_regional_combined, lcoe_results, pjs; scale_filter="Large")
save("$outputpath/fig-regional_combined_large-$opt_scaling.pdf", fig_regional_combined_large);

@info("Regional comparison plots saved (all scales + Large-specific)")

##### Learning curve plots #####
# Requires learning scenario files from smr-mcs-learning.jl
# These plots are only generated if the learning CSV files exist

# Define learning scenarios to plot (matching smr-mcs-learning.jl)
# Note: baseline is plotted as dashed reference line, not as a data point
# Three scenarios: Conservative (LR=5%), Base (LR=10%), Optimistic (LR=15%)
# κ=1.0 (FOAK-anchored), floor=nothing (unlimited learning)

learning_scenarios_conservative = [
    (1, "LR05_N1_k100"),
    (2, "LR05_N2_k100"),
    (4, "LR05_N4_k100"),
    (6, "LR05_N6_k100"),
    (12, "LR05_N12_k100")
]

learning_scenarios_base = [
    (1, "LR10_N1_k100"),
    (2, "LR10_N2_k100"),
    (4, "LR10_N4_k100"),
    (6, "LR10_N6_k100"),
    (12, "LR10_N12_k100")
]

learning_scenarios_optimistic = [
    (1, "LR15_N1_k100"),
    (2, "LR15_N2_k100"),
    (4, "LR15_N4_k100"),
    (6, "LR15_N6_k100"),
    (12, "LR15_N12_k100")
]

# Check if at least the baseline learning file exists
baseline_file = "$outputpath/mcs-lcoe_summary-$opt_scaling-baseline.csv"
if isfile(baseline_file)
    @info("Learning scenario files found - generating learning curve plots")

    # Conservative scenario (LR=5%)
    @info("Generating conservative learning curves (LR=5%)")
    fig_learning_comparison_conservative = learning_curve_comparison_plot(outputpath, opt_scaling, learning_scenarios_conservative, pjs)
    save("$outputpath/fig-learning_curve_by_scale-$opt_scaling-LR05.pdf", fig_learning_comparison_conservative)

    fig_learning_overall_conservative = learning_curve_plot(outputpath, opt_scaling, learning_scenarios_conservative)
    save("$outputpath/fig-learning_curve_overall-$opt_scaling-LR05.pdf", fig_learning_overall_conservative)

    # Base scenario (LR=10%)
    @info("Generating base learning curves (LR=10%)")
    fig_learning_comparison_base = learning_curve_comparison_plot(outputpath, opt_scaling, learning_scenarios_base, pjs)
    save("$outputpath/fig-learning_curve_by_scale-$opt_scaling-LR10.pdf", fig_learning_comparison_base)

    fig_learning_overall_base = learning_curve_plot(outputpath, opt_scaling, learning_scenarios_base)
    save("$outputpath/fig-learning_curve_overall-$opt_scaling-LR10.pdf", fig_learning_overall_base)

    # Optimistic scenario (LR=15%)
    @info("Generating optimistic learning curves (LR=15%)")
    fig_learning_comparison_optimistic = learning_curve_comparison_plot(outputpath, opt_scaling, learning_scenarios_optimistic, pjs)
    save("$outputpath/fig-learning_curve_by_scale-$opt_scaling-LR15.pdf", fig_learning_comparison_optimistic)

    fig_learning_overall_optimistic = learning_curve_plot(outputpath, opt_scaling, learning_scenarios_optimistic)
    save("$outputpath/fig-learning_curve_overall-$opt_scaling-LR15.pdf", fig_learning_overall_optimistic)

    @info("Learning curve plots saved for all three scenarios (conservative, base, optimistic)")
else
    @info("Learning scenario files not found - skipping learning curve plots")
    @info("Run smr-mcs-learning.jl first to generate learning scenario data")
end

##### WACC sensitivity analysis #####
# Shows how median LCOE varies with discount rate for different reactor scales

@info("Generating WACC sensitivity plot")

# Define WACC range for sensitivity analysis (3% to 12% in 1% steps)
wacc_sensitivity_range = 0.03:0.01:0.12

# Generate the plot (uses fewer simulations for speed: 5000 per WACC point)
fig_wacc_sensitivity = Base.invokelatest(
    wacc_sensitivity_plot,
    pjs,
    opt_scaling,
    wacc_sensitivity_range,
    electricity_price_mean,
    construction_time_ranges;
    n=5000  # Reduced from 10000 for faster computation
)

save("$outputpath/fig-wacc_sensitivity-$opt_scaling.pdf", fig_wacc_sensitivity)
@info("WACC sensitivity plot saved")
