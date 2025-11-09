##### plots #####
# Standalone plotting script - loads results from CSV files
# Compatible with split workflow (run_1_mcs.jl, run_2_sensitivity.jl, run_3_shapley.jl)

using Pkg
Pkg.activate(pwd())
using CairoMakie, CSV, DataFrames, Statistics

inputpath = "_input"
outputpath = "_output"

include("functions.jl")
include("functions_plots.jl")
include("data.jl")

##### Load results from CSV files #####

# Configuration: Select scaling method
opts_scaling = ["manufacturer", "roulstone", "rothwell", "uniform", "carelli"]
local_scaling_index = 3  # 3 = rothwell

if @isdefined(par_job) == true
    opt_scaling = opts_scaling[par_job]
    @info("Cluster job mode: using scaling option $opt_scaling")
else
    opt_scaling = opts_scaling[local_scaling_index]
    @info("Interactive mode: using scaling option $opt_scaling")
end

@info("Loading results from _output/ directory for scaling: $opt_scaling")

# Load MCS results
npv_results = CSV.read("$outputpath/mcs-npv_results-$opt_scaling.csv", DataFrame)
lcoe_results = CSV.read("$outputpath/mcs-lcoe_results-$opt_scaling.csv", DataFrame)

# Declare variables outside try blocks for proper scoping
si_npv_results = nothing
si_lcoe_results = nothing
shapley_npv_results = nothing
shapley_lcoe_results = nothing

# Load Sensitivity Index results (if they exist)
try
    global si_npv_results = CSV.read("$outputpath/si-npv_results-$opt_scaling.csv", DataFrame)
    global si_lcoe_results = CSV.read("$outputpath/si-lcoe_results-$opt_scaling.csv", DataFrame)
    @info("Loaded sensitivity index results")
catch e
    @warn("Could not load sensitivity index results: $e")
    @warn("Run run_2_sensitivity.jl first if you need SI plots")
end

# Load Shapley results (if they exist)
try
    global shapley_npv_results = CSV.read("$outputpath/shapley-npv_results-$opt_scaling.csv", DataFrame)
    global shapley_lcoe_results = CSV.read("$outputpath/shapley-lcoe_results-$opt_scaling.csv", DataFrame)
    @info("Loaded Shapley results")
catch e
    @warn("Could not load Shapley results: $e")
    @warn("Run run_3_shapley.jl first if you need Shapley plots")
end

@info("Data loaded successfully. Generating plots...")
@info("="^80)

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
carelli = hlines!(ax_theory, [0.20], label = "Carelli", linewidth = 3, color = :orange, linestyle = :dash);
hlines!(ax_theory, [0], color = :gray);

Legend(fig_theory[1, 1],
    [roulstone, rothwell, carelli],
    [L"β^\text{Roulstone}", L"β^\text{Rothwell}", L"β^\text{Carelli}"],
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
# OPTIONAL: These plots compare Roulstone vs Rothwell scaling methods
# WARNING: Requires regenerating investment data (300k additional simulations)
# Uncomment only if you need methodological comparison between scaling approaches

# hist_invest = Figure();
#
# for i in 1:3, j in 1:5
#     hist_invest_plot(n, wacc, electricity_price_mean, pjs[j+5*(i-1)], i, j, hist_invest)
# end
#
# Legend(hist_invest[4,1:5],
#     [roulstone, rothwell],
#     ["Roulstone", "Rothwell"],
#     framevisible = false, orientation = :horizontal)
#
# hist_invest
# save("$outputpath/fig-histogram_investment.pdf", hist_invest);

@info("Skipping histogram plots (commented out - saves 300k simulations). Uncomment if needed for scaling method comparison.")

##### probability density plots for comparison of estimation approaches #####
# OPTIONAL: Same as histograms above - only needed for methodological comparison
# WARNING: Requires regenerating investment data (300k additional simulations)

# density_invest = Figure();
#
# for i in 1:3, j in 1:5
#     density_invest_plot(n, wacc, electricity_price_mean, pjs[j+5*(i-1)], i, j, density_invest)
# end
#
# Legend(density_invest[4,1:5],
#     [roulstone, rothwell],
#     ["Roulstone", "Rothwell"],
#     framevisible = false, orientation = :horizontal)
#
# density_invest
# save("$outputpath/fig-density_investment.pdf", density_invest);

@info("Skipping density plots (commented out - saves 300k simulations). Uncomment if needed for scaling method comparison.");

##### boxplots Monte Carlo simulation results #####
# requires results for all 15 reactor concepts

fig_mcs_npv = mcs_plot(npv_results, "NPV", "[USD/MW]")
fig_mcs_lcoe = mcs_plot(lcoe_results, "LCOE", "[USD/MWh]")

save("$outputpath/fig-mcs_npv-$opt_scaling.pdf", fig_mcs_npv);
save("$outputpath/fig-mcs_lcoe-$opt_scaling.pdf", fig_mcs_lcoe);

##### heatmaps sensitivity indices #####
# requires sensitvity results
# Note: si_plot() was replaced by si_plot_by_scale() below for better organization

if !isnothing(si_npv_results) && !isnothing(si_lcoe_results)
    @info("Generating sensitivity index plots")
    # New: Grouped by scale (Micro/SMR/Large)
    # Use Base.invokelatest to avoid Julia 1.12 world age issues
    fig_si_npv_by_scale = Base.invokelatest(si_plot_by_scale, si_npv_results, "NPV Sensitivity Indices", pjs)
    fig_si_lcoe_by_scale = Base.invokelatest(si_plot_by_scale, si_lcoe_results, "LCOE Sensitivity Indices", pjs)

    save("$outputpath/fig-si_npv_by_scale-$opt_scaling.pdf", fig_si_npv_by_scale);
    save("$outputpath/fig-si_lcoe_by_scale-$opt_scaling.pdf", fig_si_lcoe_by_scale);
    @info("  - fig-si_npv_by_scale-$opt_scaling.pdf")
    @info("  - fig-si_lcoe_by_scale-$opt_scaling.pdf")
else
    @warn("Skipping sensitivity index plots (data not available)")
end

##### heatmaps Shapley effects #####
# Requires Shapley results from run_simulation.jl
# Shapley effects properly handle correlated inputs (WACC × Construction Time)

if !isnothing(shapley_npv_results) && !isnothing(shapley_lcoe_results)
    @info("Generating Shapley sensitivity plots")

    # Create heatmap plots grouped by scale
    fig_shapley_npv_by_scale = Base.invokelatest(shapley_plot_by_scale, shapley_npv_results, "NPV Shapley Effects", pjs)
    fig_shapley_lcoe_by_scale = Base.invokelatest(shapley_plot_by_scale, shapley_lcoe_results, "LCOE Shapley Effects", pjs)

    save("$outputpath/fig-shapley_npv_by_scale-$opt_scaling.pdf", fig_shapley_npv_by_scale);
    save("$outputpath/fig-shapley_lcoe_by_scale-$opt_scaling.pdf", fig_shapley_lcoe_by_scale);

    @info("Shapley sensitivity plots saved")
    @info("  - fig-shapley_npv_by_scale-$opt_scaling.pdf")
    @info("  - fig-shapley_lcoe_by_scale-$opt_scaling.pdf")
else
    @warn("Skipping Shapley plots (data not available)")
end

##### lcoe comparison plot #####
# Original structure preserved - adding Large and Micro to SMR results

# Load LCOE summary statistics
lcoe_summary = CSV.read("$outputpath/mcs-lcoe_summary-$opt_scaling.csv", DataFrame)

# Read external LCOE data (renewables, conventionals)
lcoe_dat = CSV.read("$inputpath/lcoe_data.csv", DataFrame)

# Calculate LCOE ranges for simulation results (Large/SMR/Micro by type)
sim_lcoe_data = DataFrame(technology=String[], lower_bound=Float64[], upper_bound=Float64[])
reactor_types = String[]  # Store reactor types separately for color mapping

# Reactor groups in display order (BOTTOM to TOP for plotting)
# Plot shows: Renewables → Conventionals → SMR → Micro → Large (bottom to top)
reactor_groups = [
    # SMR section (OLD structure - 3 rows with BWR & PWR combined)
    ("SMR", "SFR", "SFR SMRs"),
    ("SMR", "HTR", "HTR SMRs"),
    ("SMR", "BWR+PWR", "BWR & PWR SMRs"),  # Special: combines both types

    # Micro section (NEW - 3 rows)
    ("Micro", "PWR", "Micro PWR"),
    ("Micro", "HTR", "Micro HTR"),
    ("Micro", "SFR", "Micro SFR"),

    # Large section (NEW - 3 rows)
    ("Large", "PWR+BWR", "Large PWR & BWR"),  # Special: combines both types
    ("Large", "HTR", "Large HTR"),
    ("Large", "SFR", "Large SFR")
]

for (scale, rtype, label) in reactor_groups
    matching_reactors = String[]
    for pj in pjs
        # Handle special combined rows
        if rtype == "BWR+PWR" && pj.scale == scale && (pj.type == "BWR" || pj.type == "PWR")
            push!(matching_reactors, pj.name)
        elseif rtype == "PWR+BWR" && pj.scale == scale && (pj.type == "PWR" || pj.type == "BWR")
            push!(matching_reactors, pj.name)
        elseif pj.scale == scale && pj.type == rtype
            push!(matching_reactors, pj.name)
        end
    end

    if !isempty(matching_reactors)
        reactor_lcoes = []
        for reactor in matching_reactors
            # Find reactor in the variable column
            idx = findfirst(lcoe_summary.variable .== reactor)
            if !isnothing(idx)
                push!(reactor_lcoes, (lcoe_summary.q25[idx], lcoe_summary.q75[idx]))
            end
        end
        if !isempty(reactor_lcoes)
            lower = minimum([x[1] for x in reactor_lcoes])
            upper = maximum([x[2] for x in reactor_lcoes])
            push!(sim_lcoe_data, (label, lower, upper))
            # Store reactor type separately for color assignment
            reactor_type = rtype == "BWR+PWR" || rtype == "PWR+BWR" ? "PWR" : rtype
            push!(reactor_types, reactor_type)
        end
    end
end

# Combine external data with simulation results
lcoe_plot_data = vcat(
    select(lcoe_dat, [:technology, :lower_bound, :upper_bound]),
    sim_lcoe_data
)

# Define plot styling
plot_scaling = if opt_scaling == "manufacturer"
    "Manufacturer"
elseif opt_scaling == "roulstone"
    "Roulstone"
elseif opt_scaling == "rothwell"
    "Rothwell"
elseif opt_scaling == "uniform"
    "Uniform"
elseif opt_scaling == "carelli"
    "Carelli"
else
    opt_scaling
end

xlabel = "[USD/MWh]"
yticks = lcoe_plot_data[!,:technology]

# Color scheme with reactor type variation
# Renewables (rows 1-8): purple
# Conventionals (rows 9-12): dark blue
# Nuclear reactors: color by type (SFR=teal, HTR=yellow, PWR/BWR=orangered)
n_external = nrow(lcoe_dat)
n_renewables = 8
n_conventionals = n_external - n_renewables  # Should be 4

# Define color mapping by reactor type
function get_reactor_color(reactor_type)
    if reactor_type == "SFR"
        return :teal
    elseif reactor_type == "HTR"
        return :gold
    elseif reactor_type == "PWR"
        return :orangered
    else
        return :gray  # Fallback
    end
end

# Build color vector
nuclear_colors = [get_reactor_color(rt) for rt in reactor_types]

col = vcat(
    fill(:purple, n_renewables),       # Renewables (purple)
    fill(:darkblue, n_conventionals),  # Conventionals (dark blue)
    nuclear_colors                      # Nuclear (color by type)
)

fig_lcoe_comparison = Figure()
ax_lcoe = Axis(fig_lcoe_comparison[1,1],
               yticks = (1:length(yticks), yticks),
               xscale = log10,
               xlabel = xlabel)

xlims!(10, 25000)

rangebars!(ax_lcoe, 1:length(yticks), lcoe_plot_data[!,:lower_bound], lcoe_plot_data[!,:upper_bound],
           linewidth = 6, whiskerwidth = 8, direction = :x, color = col)

# Dividing lines (4 red dashed lines)
# Line 1: Between renewables and conventionals (after row 8)
renewables_end = 8
# Line 2: Between conventionals and SMR (after row 12 = all LAZARD data)
conventionals_end = n_external  # Should be 12
# Line 3: Between SMR and Micro (after 3 SMR rows)
smr_end = n_external + 3  # Row 15
# Line 4: Between Micro and Large (after 3 Micro rows)
micro_end = n_external + 6  # Row 18

hlines!(ax_lcoe, [renewables_end + 0.5, conventionals_end + 0.5, smr_end + 0.5, micro_end + 0.5],
        linestyle = :dash, color = :red, linewidth = 1.5)

# Section labels (rotated 90°, positioned at center of each section)
renewables_center = renewables_end / 2  # Center of rows 1-8 = row 4
conventionals_center = renewables_end + (conventionals_end - renewables_end) / 2  # Center of rows 9-12 = row 10.5
smr_center = n_external + 1.5  # Center of 3 SMR rows (rows 13-15) = row 14
micro_center = smr_end + 1.5  # Center of 3 Micro rows (rows 16-18) = row 17
large_center = micro_end + 1.5  # Center of 3 Large rows (rows 19-21) = row 20

text!([15000, 15000, 15000, 15000, 15000],
      [renewables_center, conventionals_center, smr_center, micro_center, large_center],
      text = ["Renewables\n(LAZARD)",
              "Conventionals\n(LAZARD)",
              "SMR\n($plot_scaling)",
              "Micro\n($plot_scaling)",
              "Large\n($plot_scaling)"],
      align = (:center, :center),
      justification = :center,
      rotation = π/2,
      fontsize = 10)

# Add values
text!(lcoe_plot_data[!,:lower_bound], 1:length(yticks),
      text = string.(round.(Int, lcoe_plot_data[!,:lower_bound])),
      align = (:right, :center), offset = (-10,0), fontsize = 9)
text!(lcoe_plot_data[!,:upper_bound], 1:length(yticks),
      text = string.(round.(Int, lcoe_plot_data[!,:upper_bound])),
      align = (:left, :center), offset = (10,0), fontsize = 9)

Label(fig_lcoe_comparison[1, 1, Top()], "LCOE Comparison",
      font = "Noto Sans Bold", padding = (0, 6, 6, 0))

fig_lcoe_comparison
save("$outputpath/fig-lcoe_comparison-$opt_scaling.pdf", fig_lcoe_comparison)
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

# Define WACC bin centers for sensitivity analysis
# Main simulation uses WACC range [4%, 10%], so we bin within this range
# Bins: 4%, 5%, 6%, 7%, 8%, 9%, 10% (each bin is ±0.5% around center)
wacc_bin_centers = 4:1:10

# NEW APPROACH: Bin existing results instead of re-running simulations
# This eliminates 720k simulation re-runs (6× speedup)
fig_wacc_sensitivity = Base.invokelatest(
    wacc_sensitivity_plot,
    outputpath,
    opt_scaling,
    pjs_dat,  # DataFrame with reactor metadata
    wacc_bin_centers
)

save("$outputpath/fig-wacc_sensitivity-$opt_scaling.pdf", fig_wacc_sensitivity)
@info("WACC sensitivity plot saved (0 new simulations, used existing data)")

##### OCC vs Year plot (Large reactors only) #####
# Shows Western vs Asian reactors with linear trend lines

@info("Generating OCC vs Year plot for Large reactors")

# Check if year column exists
if :year in propertynames(pjs_dat)
    # Filter to Large reactors only
    df_large = filter(row -> row.scale == "Large", pjs_dat)

    if nrow(df_large) > 0
        # Calculate OCC (USD/kW) from investment (USD/MW)
        df_large.occ = df_large.investment ./ 1000

        # Group by Western vs Asian
        western_regions = ["Western / Developed", "Eastern Europe", "South America"]
        asian_regions = ["Emerging Asia"]

        df_western = filter(row -> row.region in western_regions, df_large)
        df_asian = filter(row -> row.region in asian_regions, df_large)

        # Create plot
        fig_occ_year = Figure(size=(900, 600))
        ax_occ = Axis(fig_occ_year[1, 1],
                     xlabel = "Grid Connection Year",
                     ylabel = "OCC (USD2020/kW)",
                     title = "Overnight Capital Cost vs Year (Large Reactors)")

        # Plot Western reactors (blue)
        if nrow(df_western) > 0
            scatter!(ax_occ, df_western.year, df_western.occ,
                    markersize = 12, color = :blue, strokewidth = 1, strokecolor = :black,
                    label = "Western")

            # Add linear trend line (simple linear regression)
            if nrow(df_western) > 1
                x_w = Float64.(df_western.year)
                y_w = Float64.(df_western.occ)
                n_w = length(x_w)
                slope_w = (n_w * sum(x_w .* y_w) - sum(x_w) * sum(y_w)) / (n_w * sum(x_w.^2) - sum(x_w)^2)
                intercept_w = mean(y_w) - slope_w * mean(x_w)

                year_range_w = minimum(x_w):maximum(x_w)
                occ_pred_w = intercept_w .+ slope_w .* year_range_w

                lines!(ax_occ, year_range_w, occ_pred_w,
                      color = :blue, linewidth = 2, linestyle = :dash)
            end
        end

        # Plot Asian reactors (orange)
        if nrow(df_asian) > 0
            scatter!(ax_occ, df_asian.year, df_asian.occ,
                    markersize = 12, color = :orange, strokewidth = 1, strokecolor = :black,
                    label = "Asian")

            # Add linear trend line
            if nrow(df_asian) > 1
                x_a = Float64.(df_asian.year)
                y_a = Float64.(df_asian.occ)
                n_a = length(x_a)
                slope_a = (n_a * sum(x_a .* y_a) - sum(x_a) * sum(y_a)) / (n_a * sum(x_a.^2) - sum(x_a)^2)
                intercept_a = mean(y_a) - slope_a * mean(x_a)

                year_range_a = minimum(x_a):maximum(x_a)
                occ_pred_a = intercept_a .+ slope_a .* year_range_a

                lines!(ax_occ, year_range_a, occ_pred_a,
                      color = :orange, linewidth = 2, linestyle = :dash)
            end
        end

        axislegend(ax_occ, position = :rt)
        ax_occ.xgridvisible = true
        ax_occ.ygridvisible = true

        save("$outputpath/fig-occ_vs_year-large.pdf", fig_occ_year)
        @info("OCC vs Year plot saved")
    else
        @info("No Large reactors found - skipping OCC vs Year plot")
    end
else
    @info("Year column not found in reactor_data.csv - skipping OCC vs Year plot")
    @info("Run convert_to_csv.jl to regenerate CSV with year column")
end

##### Summary Tables #####
# Generate reactor dataset summary statistics tables

@info("Generating summary statistics tables")

# Check if PrettyTables is available for nice formatting
const PRETTYTABLES_AVAILABLE = try
    using PrettyTables
    true
catch
    false
end

if !PRETTYTABLES_AVAILABLE
    @warn "PrettyTables not installed. Tables will use basic formatting. Install with: import Pkg; Pkg.add(\"PrettyTables\")"
end

# Read reactor data (already available as pjs_dat)
total_reactors = nrow(pjs_dat)

println("\n" * "="^80)
println("REACTOR DATASET SUMMARY STATISTICS")
println("="^80)
println("Total reactors: $total_reactors")
println()

# ========== TABLE 1: BY SCALE ==========
scale_data = DataFrame(
    Category = String[],
    Range = String[],
    Count = Int[],
    Share = String[],
    Description = String[]
)

scale_order = ["Large", "SMR", "Micro"]
scale_ranges = Dict(
    "Large" => ">300 MWe",
    "SMR" => "50-300 MWe",
    "Micro" => "<50 MWe"
)
scale_descriptions = Dict(
    "Large" => "Conventional Gen III/III+ units",
    "SMR" => "Modular designs under development",
    "Micro" => "Small-scale/off-grid designs"
)

scale_counts = combine(groupby(pjs_dat, :scale), nrow => :count)
for scale in scale_order
    row = filter(r -> r.scale == scale, scale_counts)
    if nrow(row) > 0
        count = row[1, :count]
        share = "$(round(count / total_reactors * 100, digits=1))%"
        push!(scale_data, (scale, scale_ranges[scale], count, share, scale_descriptions[scale]))
    end
end

println("="^80)
println("TABLE 1: REACTORS BY SCALE")
println("="^80)
if PRETTYTABLES_AVAILABLE
    pretty_table(scale_data,
         backend=:text,
         alignment=[:l, :l, :r, :r, :l])
else
    println(scale_data)
end

# ========== TABLE 2: BY TYPE ==========
type_data = DataFrame(
    Type = String[],
    Full_Name = String[],
    Count = Int[],
    Share = String[],
    Description = String[]
)

type_full_names = Dict(
    "PWR" => "Pressurized Water Reactor",
    "BWR" => "Boiling Water Reactor",
    "SFR" => "Sodium-cooled Fast Reactor",
    "HTR" => "High-Temperature Gas-cooled Reactor"
)

type_descriptions = Dict(
    "PWR" => "Dominant type in SMR and Large",
    "BWR" => "GE, Hitachi designs",
    "SFR" => "BN, CEFR, ARC designs",
    "HTR" => "Fort St. Vrain, HTR-PM, EM2"
)

type_counts = combine(groupby(pjs_dat, :type), nrow => :count)
sort!(type_counts, :count, rev=true)

for row in eachrow(type_counts)
    rtype = row.type
    count = row.count
    share = "$(round(count / total_reactors * 100, digits=1))%"
    full_name = get(type_full_names, rtype, rtype)
    desc = get(type_descriptions, rtype, "")
    push!(type_data, (rtype, full_name, count, share, desc))
end

println("\n" * "="^80)
println("TABLE 2: REACTORS BY TYPE")
println("="^80)
if PRETTYTABLES_AVAILABLE
    pretty_table(type_data,
         backend=:text,
         alignment=[:l, :l, :r, :r, :l])
else
    println(type_data)
end

# ========== TABLE 3: BY REGION ==========
region_data = DataFrame(
    Region = String[],
    Count = Int[],
    Share = String[],
    Countries = String[],
    Description = String[]
)

region_countries = Dict(
    "Western / Developed" => "USA, UK, France, Canada, Japan, S. Korea",
    "Emerging Asia" => "China, Russia, India, Pakistan",
    "South America" => "Argentina, Brazil",
    "Middle East / Africa" => "UAE, Saudi Arabia, Egypt, S. Africa"
)

region_descriptions = Dict(
    "Western / Developed" => "OECD vendor-led projects",
    "Emerging Asia" => "State-driven programs",
    "South America" => "Regional development",
    "Middle East / Africa" => "Emerging markets"
)

region_counts = combine(groupby(pjs_dat, :region), nrow => :count)
sort!(region_counts, :count, rev=true)

for row in eachrow(region_counts)
    region = row.region
    count = row.count
    share = "$(round(count / total_reactors * 100, digits=1))%"
    countries = get(region_countries, region, "Various")
    desc = get(region_descriptions, region, "")
    push!(region_data, (region, count, share, countries, desc))
end

println("\n" * "="^80)
println("TABLE 3: REACTORS BY GEOGRAPHIC REGION")
println("="^80)
if PRETTYTABLES_AVAILABLE
   pretty_table(region_data,
         backend=:text,
         alignment=[:l, :r, :r, :l, :l])
else
    println(region_data)
end

println("\n" * "="^80)

# Save tables to CSV in _output/ directory
CSV.write("$outputpath/summary_by_scale.csv", scale_data)
CSV.write("$outputpath/summary_by_type.csv", type_data)
CSV.write("$outputpath/summary_by_region.csv", region_data)

@info("Summary tables saved to $outputpath/")
@info("  - summary_by_scale.csv")
@info("  - summary_by_type.csv")
@info("  - summary_by_region.csv")
