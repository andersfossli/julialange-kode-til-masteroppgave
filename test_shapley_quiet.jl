##### Test Shapley Quiet Mode - Run on 2 reactors only #####

using Pkg
Pkg.activate(pwd())
using Statistics

inputpath = "_input"
outputpath = "_output"

@info("="^80)
@info("TESTING SHAPLEY QUIET MODE (2 reactors only)")
@info("="^80)

##### load functions #####
@info("Loading functions")
include("functions.jl")

##### load project data #####
@info("Loading data")
include("data.jl")

@info("Total reactors available: $(length(pjs))")
@info("Testing with first 2 reactors only")

##### simulation parameters #####
n = Int64(100000)
electricity_price_mean = mean([52.2, 95.8])
wacc = [0.04, 0.10]

construction_time_ranges = Dict(
    "Micro" => [3, 7],
    "SMR"   => [3, 7],
    "Large" => [5, 12]
)

# scaling parameter
scaling = [0.4, 0.7]  # Thesis: uniform distribution

opts_scaling = ["manufacturer", "roulstone", "rothwell", "uniform", "carelli"]
local_scaling_index = 2

if @isdefined(par_job) == true
    opt_scaling = opts_scaling[par_job]
else
    opt_scaling = opts_scaling[local_scaling_index]
end

@info("Using scaling option: $opt_scaling")

##### TEST WITH 2 REACTORS ONLY #####
test_reactors = [1, 2]  # First two reactors

for p in test_reactors
    @info("="^80)
    @info("Testing reactor $(p)/2: $(pjs[p].name)")
    @info("="^80)

    construction_time_range = construction_time_ranges[pjs[p].scale]

    # Run with quiet=true
    @info("Running Shapley analysis with QUIET MODE enabled...")
    shapley_results = shapley_sensitivity_index(
        opt_scaling, n, wacc, electricity_price_mean, pjs[p];
        construction_time_range=construction_time_range,
        quiet=true
    )

    @info("Results for $(pjs[p].name):")
    @info("  Sh_NPV:  $(shapley_results.sh_npv)")
    @info("  Sh_LCOE: $(shapley_results.sh_lcoe)")
    @info("Completed $(p)/2 reactors")

    GC.gc()
end

@info("="^80)
@info("TEST COMPLETE - Check output above!")
@info("="^80)
@info("If you see minimal logging (no coalition details), quiet mode is working!")
