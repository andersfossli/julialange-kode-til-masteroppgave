##### plots #####
# Standalone plotting script - reads results from CSV files
# Can be run independently after simulations have been completed

##### Load dependencies #####
using Pkg
Pkg.activate(pwd())
using CairoMakie
using CSV
using DataFrames
using Statistics

inputpath = "_input"
outputpath = "_output"

@info("="^80)
@info("GENERATING PLOTS FROM SIMULATION RESULTS")
@info("="^80)

##### Load core functions and data #####
@info("Loading functions and data")
include("functions.jl")
include("data.jl")
include("config.jl")
include("functions_plots.jl")

##### Load simulation results from CSV files #####
@info("Loading simulation results from CSV files")

# Check if required files exist
required_files = [
    "$outputpath/mcs-npv_results-$opt_scaling.csv",
    "$outputpath/mcs-npv_summary-$opt_scaling.csv",
    "$outputpath/mcs-lcoe_results-$opt_scaling.csv",
    "$outputpath/mcs-lcoe_summary-$opt_scaling.csv"
]

missing_files = filter(f -> !isfile(f), required_files)
if !isempty(missing_files)
    @error("Missing required result files. Please run simulations first:")
    for f in missing_files
        @error("  - $f")
    end
    error("Cannot generate plots without simulation results")
end

# Load Monte Carlo results
npv_results = CSV.read("$outputpath/mcs-npv_results-$opt_scaling.csv", DataFrame)
npv_summary = CSV.read("$outputpath/mcs-npv_summary-$opt_scaling.csv", DataFrame)
lcoe_results = CSV.read("$outputpath/mcs-lcoe_results-$opt_scaling.csv", DataFrame)
lcoe_summary = CSV.read("$outputpath/mcs-lcoe_summary-$opt_scaling.csv", DataFrame)

@info("✓ Loaded Monte Carlo results")

# Load sensitivity indices (optional - some plots will be skipped if missing)
if isfile("$outputpath/si-npv_results-$opt_scaling.csv")
    si_npv_results = CSV.read("$outputpath/si-npv_results-$opt_scaling.csv", DataFrame)
    si_lcoe_results = CSV.read("$outputpath/si-lcoe_results-$opt_scaling.csv", DataFrame)
    @info("✓ Loaded Sobol sensitivity results")
else
    @warn("Sobol sensitivity results not found - skipping sensitivity plots")
    @warn("Run run_2_sensitivity.jl to generate these results")
    si_npv_results = nothing
    si_lcoe_results = nothing
end

# Load Shapley results (optional)
if isfile("$outputpath/shapley-npv_results-$opt_scaling.csv")
    shapley_npv_results = CSV.read("$outputpath/shapley-npv_results-$opt_scaling.csv", DataFrame)
    shapley_lcoe_results = CSV.read("$outputpath/shapley-lcoe_results-$opt_scaling.csv", DataFrame)
    @info("✓ Loaded Shapley sensitivity results")
else
    @warn("Shapley sensitivity results not found - some plots may be skipped")
    @warn("Run run_3_shapley.jl to generate these results")
    shapley_npv_results = nothing
    shapley_lcoe_results = nothing
end

@info("="^80)
@info("STARTING PLOT GENERATION")
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

# NPV plots disabled - not needed for thesis
# fig_mcs_npv = mcs_plot(npv_results, "NPV", "[EUR2025/MW]", pjs)
# save("$outputpath/fig-mcs_npv-$opt_scaling.pdf", fig_mcs_npv);

fig_mcs_lcoe = mcs_plot(lcoe_results, "LCOE", "[EUR2025/MWh]", pjs)
save("$outputpath/fig-mcs_lcoe-$opt_scaling.pdf", fig_mcs_lcoe);

##### heatmaps sensitivity indices #####
# requires sensitvity results
# Note: si_plot() was replaced by si_plot_by_scale() below for better organization

if !isnothing(si_npv_results) && !isnothing(si_lcoe_results)
    @info("Generating sensitivity index plots")
    # New: Grouped by scale (Micro/SMR/Large)
    # Use Base.invokelatest to avoid Julia 1.12 world age issues

    # NPV sensitivity plots disabled - not needed for thesis
    # fig_si_npv_by_scale = Base.invokelatest(si_plot_by_scale, si_npv_results, "NPV Sensitivity Indices", pjs)
    # save("$outputpath/fig-si_npv_by_scale-$opt_scaling.pdf", fig_si_npv_by_scale);

    fig_si_lcoe_by_scale = Base.invokelatest(si_plot_by_scale, si_lcoe_results, "LCOE Sensitivity Indices", pjs)
    save("$outputpath/fig-si_lcoe_by_scale-$opt_scaling.pdf", fig_si_lcoe_by_scale);
    @info("✓ Sensitivity index plots saved")
else
    @warn("Skipping sensitivity index plots (data not available)")
end

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
xlabel = "[EUR2025/MWh]";
yticks = lcoe_plot_data[!,:technology];

# Dynamic color assignment based on actual data size
n_lcoe_data = nrow(lcoe_dat)  # Number of rows from lcoe_data.csv
n_smr_summary = 3              # BWR & PWR SMRs, HTR SMRs, SFR SMRs

# Assign colors: first n_lcoe_data rows get incrementing colors, last 3 get distinct colors
col = vcat(collect(1:n_lcoe_data), (n_lcoe_data+1):(n_lcoe_data+n_smr_summary))
colormap = [:darkgreen, :darkblue];

fig_lcoe_comparison = Figure();
ax_lcoe = Axis(fig_lcoe_comparison[1,1], yticks = (1:length(yticks), yticks), xscale = log10, xlabel = xlabel);

xlims!(10, 25000)

rangebars!(ax_lcoe, 1:length(yticks), lcoe_plot_data[!,2], lcoe_plot_data[!,3], linewidth = 6, whiskerwidth = 8, direction = :x, color = col);

# Calculate separator positions based on actual data
n_renewables = 7  # PV×3, Geothermal, Wind×3
n_conventionals = 4  # Gas-Peaking, Nuclear, Coal, Gas-CombinedCycle
renewables_end = n_renewables
conventionals_end = n_renewables + n_conventionals

hlines!(ax_lcoe, [renewables_end + 0.5, conventionals_end + 0.5], linestyle = :dash, color = :red);

# Calculate text positions for category labels
renewables_center = renewables_end / 2
conventionals_center = renewables_end + (n_conventionals / 2) + 0.5
smr_center = conventionals_end + 2

text!([15000, 15000, 16], [renewables_center, conventionals_center, smr_center];
      text = ["Renewables\n(LAZARD 2025)", "Conventionals\n(LAZARD 2025)", "SMR Tech.\n($plot_scaling)"],
      align = (:center, :center), justification = :center, rotation = π/2);

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

# Define thresholds (0 to 500 EUR/MWh in steps of 20)
lcoe_thresholds = collect(0.0:20.0:500.0)  # Float64 values required
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

# Regional combined plots disabled - not needed for thesis
# fig_regional_combined_all = Base.invokelatest(mcs_plot_regional_combined, lcoe_results, pjs; scale_filter="All")
# save("$outputpath/fig-regional_combined_all-$opt_scaling.pdf", fig_regional_combined_all);

# LARGE ONLY: Specific comparison for deployment learning analysis
@info("Generating regional plots for Large reactors only")
fig_regional_large = Base.invokelatest(mcs_plot_regional, lcoe_results, pjs; scale_filter="Large")
save("$outputpath/fig-regional_comparison_large-$opt_scaling.pdf", fig_regional_large);

# Regional combined plots disabled - not needed for thesis
# fig_regional_combined_large = Base.invokelatest(mcs_plot_regional_combined, lcoe_results, pjs; scale_filter="Large")
# save("$outputpath/fig-regional_combined_large-$opt_scaling.pdf", fig_regional_combined_large);

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

    # Learning curve "by_scale" and "overall" plots disabled - not needed for thesis
    # These were replaced by the individual reactor learning curves in the section below

    # # Conservative scenario (LR=5%)
    # @info("Generating conservative learning curves (LR=5%)")
    # fig_learning_comparison_conservative = learning_curve_comparison_plot(outputpath, opt_scaling, learning_scenarios_conservative, pjs)
    # save("$outputpath/fig-learning_curve_by_scale-$opt_scaling-LR05.pdf", fig_learning_comparison_conservative)

    # fig_learning_overall_conservative = learning_curve_plot(outputpath, opt_scaling, learning_scenarios_conservative)
    # save("$outputpath/fig-learning_curve_overall-$opt_scaling-LR05.pdf", fig_learning_overall_conservative)

    # # Base scenario (LR=10%)
    # @info("Generating base learning curves (LR=10%)")
    # fig_learning_comparison_base = learning_curve_comparison_plot(outputpath, opt_scaling, learning_scenarios_base, pjs)
    # save("$outputpath/fig-learning_curve_by_scale-$opt_scaling-LR10.pdf", fig_learning_comparison_base)

    # fig_learning_overall_base = learning_curve_plot(outputpath, opt_scaling, learning_scenarios_base)
    # save("$outputpath/fig-learning_curve_overall-$opt_scaling-LR10.pdf", fig_learning_overall_base)

    # # Optimistic scenario (LR=15%)
    # @info("Generating optimistic learning curves (LR=15%)")
    # fig_learning_comparison_optimistic = learning_curve_comparison_plot(outputpath, opt_scaling, learning_scenarios_optimistic, pjs)
    # save("$outputpath/fig-learning_curve_by_scale-$opt_scaling-LR15.pdf", fig_learning_comparison_optimistic)

    # fig_learning_overall_optimistic = learning_curve_plot(outputpath, opt_scaling, learning_scenarios_optimistic)
    # save("$outputpath/fig-learning_curve_overall-$opt_scaling-LR15.pdf", fig_learning_overall_optimistic)

    @info("Learning curve by_scale and overall plots disabled (not needed for thesis)")
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
        # Calculate OCC (EUR/kW) from investment (EUR/MW)
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
                     ylabel = "OCC (EUR2025/kW)",
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

##### NEW THESIS PLOTS #####
@info("")
@info("="^80)
@info("GENERATING ADDITIONAL THESIS PLOTS")
@info("="^80)

# Define simulation parameters for plotting functions
wacc = [0.04, 0.10]  # WACC range used in simulations
electricity_price_mean = mean([52.2, 95.8])  # Mean electricity price

# PLOT 1: LCOE Comparison Horizontal Bar Chart
@info("Generating LCOE comparison horizontal bar chart...")
fig_lcoe_horizontal = lcoe_comparison_horizontal(lcoe_results, pjs, opt_scaling)
save("$outputpath/fig-lcoe_comparison_horizontal-$opt_scaling.pdf", fig_lcoe_horizontal)
@info("✓ Saved: fig-lcoe_comparison_horizontal-$opt_scaling.pdf")

# PLOT 2: Shapley Sensitivity Heatmap (matching Sobol style)
if !isnothing(shapley_npv_results) && !isnothing(shapley_lcoe_results)
    @info("Generating Shapley sensitivity heatmap (matching Sobol style)...")
    fig_shapley_heatmap = shapley_plot_by_scale(shapley_lcoe_results, "Shapley Sensitivity Indices (LCOE)", pjs)
    save("$outputpath/fig-shapley_heatmap_threepanel-$opt_scaling.pdf", fig_shapley_heatmap)
    @info("✓ Saved: fig-shapley_heatmap_threepanel-$opt_scaling.pdf")
else
    @warn("Skipping Shapley heatmap (Shapley results not available)")
end

# PLOT 3: Learning Curves - Separate plots for each scale
@info("Generating learning curves for all scales...")

# Create separate learning curve plots for each scale
large_reactors = [pj for pj in pjs if pj.scale == "Large"]
micro_reactors = [pj for pj in pjs if pj.scale == "Micro"]
smr_reactors_only = [pj for pj in pjs if pj.scale == "SMR"]

if !isempty(large_reactors)
    fig_learning_large = learning_curves_smr(large_reactors, wacc, electricity_price_mean, opt_scaling, outputpath)
    save("$outputpath/fig-learning_curves_large-$opt_scaling.pdf", fig_learning_large)
    @info("✓ Saved: fig-learning_curves_large-$opt_scaling.pdf")
end

if !isempty(micro_reactors)
    fig_learning_micro = learning_curves_smr(micro_reactors, wacc, electricity_price_mean, opt_scaling, outputpath)
    save("$outputpath/fig-learning_curves_micro-$opt_scaling.pdf", fig_learning_micro)
    @info("✓ Saved: fig-learning_curves_micro-$opt_scaling.pdf")
end

if !isempty(smr_reactors_only)
    fig_learning_smr = learning_curves_smr(smr_reactors_only, wacc, electricity_price_mean, opt_scaling, outputpath)
    save("$outputpath/fig-learning_curves_smr-$opt_scaling.pdf", fig_learning_smr)
    @info("✓ Saved: fig-learning_curves_smr-$opt_scaling.pdf")
end

# PLOT 4: Threshold Probability Curves (3 separate figures by scale)
@info("Generating threshold probability curves...")
threshold_figures = threshold_probability_curves(lcoe_results, pjs)
for (scale, fig) in threshold_figures
    filename = "$outputpath/fig-threshold_probability_$(lowercase(scale))-$opt_scaling.pdf"
    save(filename, fig)
    @info("✓ Saved: fig-threshold_probability_$(lowercase(scale))-$opt_scaling.pdf")
end

@info("")
@info("="^80)
@info("PLOT GENERATION COMPLETE")
@info("="^80)
@info("All plots saved to $outputpath/")
