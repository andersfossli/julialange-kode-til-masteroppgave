##### LCOE Component Breakdown Analysis #####
# This script generates stacked bar charts showing the breakdown of LCOE into:
# - OCC (Overnight Construction Cost)
# - IDC (Interest During Construction)
# - Fixed O&M
# - Variable O&M + Fuel
#
# Outputs: 3 PDF files (one per reactor scale: Micro, SMR, Large)

using Pkg
Pkg.activate(pwd())
using CairoMakie, DataFrames, CSV, Statistics

# Set paths
inputpath = "_input"
outputpath = "_output"

# Load functions and project data
include("functions.jl")
include("data.jl")
include("functions_plots.jl")
include("config.jl")  # Load centralized configuration

@info "="^80
@info "LCOE COMPONENT BREAKDOWN ANALYSIS"
@info "="^80

# Configuration
n = 10000  # Number of Monte Carlo simulations (10k for quick breakdown analysis)

# WACC and electricity price ranges (matching run_1_mcs.jl)
wacc = [0.04, 0.10]
electricity_price_mean = mean([52.2, 95.8])

# Construction time ranges by scale (matching run_1_mcs.jl)
construction_time_ranges = Dict(
    "Micro" => [3, 8],   # Thesis: 3-8 years, triangular mode 5
    "SMR"   => [3, 7],   # Thesis: 3-7 years, triangular mode 5
    "Large" => [5, 13]   # Thesis: 5-13 years, triangular mode 8
)

# Scaling parameter (thesis: uniform distribution)
scaling = [0.4, 0.7]

@info "Configuration:"
@info "  Scaling method: $opt_scaling"
@info ""

# Check if component breakdown CSV files exist from run_1_mcs.jl
component_files_exist = isfile("$outputpath/mcs-lcoe_occ-$opt_scaling.csv") &&
                        isfile("$outputpath/mcs-lcoe_idc-$opt_scaling.csv") &&
                        isfile("$outputpath/mcs-lcoe_fixed_om-$opt_scaling.csv") &&
                        isfile("$outputpath/mcs-lcoe_variable_om_fuel-$opt_scaling.csv")

all_breakdown_results = Dict{String, NamedTuple}()

if component_files_exist
    @info "Found existing component breakdown data from run_1_mcs.jl"
    @info "Loading LCOE component breakdown from CSV files..."

    # Load component breakdown CSVs
    lcoe_occ = CSV.read("$outputpath/mcs-lcoe_occ-$opt_scaling.csv", DataFrame)
    lcoe_idc = CSV.read("$outputpath/mcs-lcoe_idc-$opt_scaling.csv", DataFrame)
    lcoe_fixed_om = CSV.read("$outputpath/mcs-lcoe_fixed_om-$opt_scaling.csv", DataFrame)
    lcoe_variable_om_fuel = CSV.read("$outputpath/mcs-lcoe_variable_om_fuel-$opt_scaling.csv", DataFrame)

    # Build results dictionary
    for pj in pjs
        if pj.name in names(lcoe_occ)
            all_breakdown_results[pj.name] = (
                lcoe_occ = lcoe_occ[!, pj.name],
                lcoe_idc = lcoe_idc[!, pj.name],
                lcoe_capital = lcoe_occ[!, pj.name] .+ lcoe_idc[!, pj.name],
                lcoe_fixed_om = lcoe_fixed_om[!, pj.name],
                lcoe_variable_om_fuel = lcoe_variable_om_fuel[!, pj.name],
                lcoe = lcoe_occ[!, pj.name] .+ lcoe_idc[!, pj.name] .+ lcoe_fixed_om[!, pj.name] .+ lcoe_variable_om_fuel[!, pj.name]
            )

            # Quick summary
            median_occ = median(all_breakdown_results[pj.name].lcoe_occ)
            median_idc = median(all_breakdown_results[pj.name].lcoe_idc)
            median_fixed_om = median(all_breakdown_results[pj.name].lcoe_fixed_om)
            median_variable = median(all_breakdown_results[pj.name].lcoe_variable_om_fuel)
            # FIXED: Use median of totals, not sum of medians
            median_total = median(all_breakdown_results[pj.name].lcoe)

            @info "  $(pj.name): $(round(median_total, digits=1)) EUR/MWh (OCC: $(round(median_occ/median_total*100, digits=1))%, IDC: $(round(median_idc/median_total*100, digits=1))%)"
        end
    end
else
    @info "Component breakdown files not found - running new simulations"
    @info "  Simulations: $n per reactor"
    @info "  WACC range: $(wacc[1]*100)% - $(wacc[2]*100)%"
    @info "  Mean electricity price: $(round(electricity_price_mean, digits=2)) EUR/MWh"
    @info ""
    @info "Running Monte Carlo simulations with LCOE component tracking..."

    for (idx, pj) in enumerate(pjs)
        @info "  Processing $(pj.name) ($(idx)/$(length(pjs))) - $(pj.scale) $(pj.type)..."

        # Get construction time range for this scale
        construction_time_range = construction_time_ranges[pj.scale]

        # Generate random variables
        rand_vars = gen_rand_vars(opt_scaling, n, wacc, electricity_price_mean, pj;
                                 construction_time_range=construction_time_range)

        # Run Monte Carlo simulation (now includes component tracking)
        disc_res = mc_run(n, pj, rand_vars)

        # Calculate LCOE with component decomposition
        res = npv_lcoe(disc_res, decompose=true)

        # Store results
        all_breakdown_results[pj.name] = res

        # Quick summary
        median_capital = median(res.lcoe_capital)
        median_fixed_om = median(res.lcoe_fixed_om)
        median_variable = median(res.lcoe_variable_om_fuel)
        median_total = median_capital + median_fixed_om + median_variable

        @info "    Median LCOE: $(round(median_total, digits=1)) EUR2025/MWh"
        @info "      Capital: $(round(median_capital, digits=1)) ($(round(median_capital/median_total*100, digits=1))%)"
        @info "      Fixed O&M: $(round(median_fixed_om, digits=1)) ($(round(median_fixed_om/median_total*100, digits=1))%)"
        @info "      Variable O&M+Fuel: $(round(median_variable, digits=1)) ($(round(median_variable/median_total*100, digits=1))%)"
    end
end

@info ""
@info "="^80
@info "GENERATING STACKED BAR CHARTS BY SCALE"
@info "="^80

# Create stacked bar charts (one per scale)
breakdown_figures = plot_lcoe_breakdown_by_scale(pjs, all_breakdown_results, opt_scaling)

# Save figures
for (scale, fig) in breakdown_figures
    filename = "$outputpath/fig-lcoe_breakdown_$(lowercase(scale))-$opt_scaling.pdf"
    save(filename, fig)
    @info "✓ Saved: $filename"
end

@info ""
@info "="^80
@info "LCOE BREAKDOWN ANALYSIS COMPLETE"
@info "="^80
@info "Generated $(length(breakdown_figures)) stacked bar chart(s)"
@info "All plots saved to $outputpath/"
