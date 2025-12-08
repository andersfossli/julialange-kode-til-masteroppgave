##### LCOE Breakdown - Standalone Donut Chart #####
# Creates a refined standalone donut chart with vibrant colors and borders
# showing LCOE component breakdown for aggregated SMR reactors.
#
# Outputs: Single PDF with polished donut visualization

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
@info "LCOE BREAKDOWN - STANDALONE DONUT CHART"
@info "="^80

# Analysis for SMR category
scale = "SMR"
wacc = [0.04, 0.10]  # WACC range (matching run_1_mcs.jl)
electricity_price_mean = mean([52.2, 95.8])
construction_time_range = [3, 8]  # SMR construction time range
scaling = [0.4, 0.7]  # Scaling parameter range

@info "Analysis Configuration:"
@info "  Reactor scale: $scale (aggregated across all $(scale) reactors)"
@info "  WACC range: $(wacc[1]*100)% - $(wacc[2]*100)%"
@info "  Mean electricity price: $(round(electricity_price_mean, digits=2)) EUR/MWh"
@info "  Construction time range: $(construction_time_range[1])-$(construction_time_range[2]) years"
@info ""

# Generate standalone donut chart
fig = plot_lcoe_donut_standalone(pjs, scale, wacc, electricity_price_mean, opt_scaling, construction_time_range)

filename = "$outputpath/fig-lcoe_donut_$(scale)_standalone.pdf"
save(filename, fig)

@info ""
@info "="^80
@info "LCOE DONUT CHART COMPLETE"
@info "="^80
@info "✓ Saved: $filename"
@info ""
@info "Standalone donut chart with vibrant colors and borders."
