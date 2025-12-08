##### IDC Sensitivity Analysis - AGGREGATED #####
# Shows impact of construction delays on capital costs (OCC vs IDC)
# using aggregated reactor data for clearer academic storytelling.
#
# This demonstrates the general principle that construction time delays
# significantly increase Interest During Construction (IDC), impacting
# total capital costs and LCOE for nuclear projects.
#
# Outputs: Single aggregated PDF per reactor scale showing the principle

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
@info "IDC SENSITIVITY ANALYSIS (AGGREGATED)"
@info "="^80

# Analysis for SMR category (most relevant for thesis)
scale = "SMR"
wacc = 0.07  # 7% - typical WACC for nuclear projects
base_time = 3    # On-time: 3 years
delayed_time = 6  # Delayed: 6 years (2× overrun)

@info "Analysis Configuration:"
@info "  Reactor scale: $scale (aggregated across all $(scale) reactors)"
@info "  WACC: $(round(wacc*100, digits=0))%"
@info "  Construction time scenarios:"
@info "    On-time: $base_time years"
@info "    Delayed: $delayed_time years (+$(delayed_time - base_time) years)"
@info ""

# Generate aggregated IDC analysis
fig = plot_idc_aggregate(pjs, scale, wacc, base_time, delayed_time, opt_scaling)

filename = "$outputpath/fig-idc_sensitivity_$(scale)_aggregated.pdf"
save(filename, fig)

@info ""
@info "="^80
@info "IDC SENSITIVITY ANALYSIS COMPLETE"
@info "="^80
@info "✓ Saved: $filename"
@info ""
@info "This aggregated plot demonstrates the general principle of how"
@info "construction delays impact nuclear economics through IDC accumulation."
