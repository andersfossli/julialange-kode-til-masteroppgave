##### Test Shapley analysis on ONE reactor to diagnose crash #####

using Pkg
Pkg.activate(pwd())
using Statistics

inputpath = "_input"
outputpath = "_output"

@info("Loading functions")
include("functions.jl")

@info("Loading data")
include("data.jl")

# Use production parameters
n = Int64(100000)
electricity_price_mean = mean([52.2, 95.8])
wacc = [0.04, 0.10]

construction_time_ranges = Dict(
    "Micro" => [3, 7],
    "SMR"   => [3, 7],
    "Large" => [5, 12]
)

# Scaling parameter (required by gen_rand_vars)
scaling = [0.20, 0.75]

opt_scaling = "rothwell"

@info("Testing Shapley analysis on FIRST reactor only")
@info("Reactor: $(pjs[1].name)")
@info("Scale: $(pjs[1].scale)")
@info("Parameters: n=$n, n_outer=50, n_inner=100")

try
    construction_time_range = construction_time_ranges[pjs[1].scale]

    @info("Starting Shapley analysis...")
    shapley_results = shapley_sensitivity_index(
        opt_scaling, n, wacc, electricity_price_mean, pjs[1];
        construction_time_range=construction_time_range
    )

    @info("SUCCESS! Shapley analysis completed for $(pjs[1].name)")
    @info("NPV Shapley effects:")
    println("  WACC: $(shapley_results.sh_npv.wacc)")
    println("  CT:   $(shapley_results.sh_npv.construction_time)")
    println("  LF:   $(shapley_results.sh_npv.loadfactor)")
    println("  INV:  $(shapley_results.sh_npv.investment)")

    @info("LCOE Shapley effects:")
    println("  WACC: $(shapley_results.sh_lcoe.wacc)")
    println("  CT:   $(shapley_results.sh_lcoe.construction_time)")
    println("  LF:   $(shapley_results.sh_lcoe.loadfactor)")
    println("  INV:  $(shapley_results.sh_lcoe.investment)")

catch e
    @error("FAILED! Shapley analysis crashed")
    @error("Error type: $(typeof(e))")
    @error("Error message: $e")
    rethrow(e)
end
