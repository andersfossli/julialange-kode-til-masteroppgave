##### LCOE Breakdown - Aggregated (Bar + Donut) #####
# Creates side-by-side comparison of stacked bar chart and donut chart
# showing LCOE component breakdown for aggregated SMR reactors.
#
# This helps visualize which representation (bar vs. donut) works best for thesis.
#
# Outputs: Single PDF with both visualizations

using Pkg
Pkg.activate(pwd())
using CairoMakie, Statistics

# Set paths
inputpath = "_input"
outputpath = "_output"

# Load functions and project data
include("functions.jl")
include("data.jl")
include("functions_plots.jl")
include("config.jl")  # Get opt_scaling parameter

@info "="^80
@info "LCOE BREAKDOWN - AGGREGATED (BAR + DONUT COMPARISON)"
@info "="^80

# Analysis for SMR category
scale = "SMR"
wacc = [0.04, 0.10]  # WACC range (matching run_1_mcs.jl)
electricity_price_mean = mean([52.2, 95.8])
construction_time_range = [3, 7]  # SMR construction time range (3-7 years, triangular mode 5)
scaling = [0.4, 0.7]  # Scaling parameter range

@info "Analysis Configuration:"
@info "  Reactor scale: $scale (aggregated across all $(scale) reactors)"
@info "  WACC range: $(wacc[1]*100)% - $(wacc[2]*100)%"
@info "  Mean electricity price: $(round(electricity_price_mean, digits=2)) EUR/MWh"
@info "  Construction time range: $(construction_time_range[1])-$(construction_time_range[2]) years"
@info ""

# Generate aggregated breakdown with both bar and donut charts
fig = plot_lcoe_breakdown_aggregate(pjs, scale, wacc, electricity_price_mean, opt_scaling, construction_time_range)

filename = "$outputpath/fig-lcoe_breakdown_$(scale)_aggregated_comparison.pdf"
save(filename, fig)

@info ""
@info "="^80
@info "LCOE BREAKDOWN AGGREGATE COMPLETE"
@info "="^80
@info "✓ Saved: $filename"
@info ""
@info "This plot shows both bar and donut representations for comparison."
@info "Use this to decide which visualization style works best for your thesis."
