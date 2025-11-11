##### Run ONLY Shapley Sensitivity Analysis #####
# Part 3 of 3-part workflow
# Runtime: ~13-14 hours with 100k samples, n_outer=50 (41 reactors × 20 min each)
# Outputs: shapley-npv_results-*.csv, shapley-lcoe_results-*.csv

using Pkg
Pkg.activate(pwd())
using Statistics

inputpath = "_input"
outputpath = "_output"

@info("="^80)
@info("STEP 3/3: SHAPLEY SENSITIVITY ANALYSIS")
@info("="^80)

##### load functions #####
@info("Loading functions")
include("functions.jl")

##### load project data #####
@info("Loading data")
include("data.jl")

##### simulation parameters #####

# number of Monte Carlo runs (NOTE: Shapley automatically caps base variance at 10k)
# The Shapley function uses n_outer×n_inner for the main computation,
# and caps base variance samples at min(n, 10000) for efficiency
n = Int64(100000)  # Passed but capped internally to avoid wasteful base samples

# wholesale electricity price [USD/MWh]
electricity_price_mean = mean([52.2, 95.8])

# weighted average cost of capital (WACC)
wacc = [0.04, 0.10]

# construction time ranges by scale [years]
construction_time_ranges = Dict(
    "Micro" => [3, 7],
    "SMR"   => [3, 7],
    "Large" => [5, 12]
)

# scaling parameter
scaling = [0.4, 0.7]  # Thesis: uniform distribution

# scaling options
opts_scaling = ["manufacturer", "roulstone", "rothwell", "uniform", "carelli"]
local_scaling_index = 2  # 2 = roulstone (thesis default)

# choose scaling option
if @isdefined(par_job) == true
    opt_scaling = opts_scaling[par_job]
    @info("Cluster job mode: using scaling option $opt_scaling")
else
    opt_scaling = opts_scaling[local_scaling_index]
    @info("Interactive mode: using scaling option $opt_scaling")
end

##### RUN SHAPLEY SENSITIVITY ANALYSIS ONLY #####

# Initialize Shapley results variables
shapley_npv_results = DataFrame()
shapley_npv_results.var = ["wacc", "construction_time", "loadfactor", "investment"]
shapley_lcoe_results = DataFrame()
shapley_lcoe_results.var = ["wacc", "construction_time", "loadfactor", "investment"]

@info("Starting Shapley sensitivity analysis for all reactors")
@info("Total reactors: $(length(pjs))")
@info("Expected runtime: ~$(round(length(pjs) * 20/60, digits=1)) hours ($(length(pjs)) reactors × ~20 min each)")
@info("Using QUIET MODE to reduce log file size (coalition details suppressed)")

# Run Shapley sensitivity analysis for all projects
for p in eachindex(pjs)
    @info("Running Shapley sensitivity for reactor $(p)/$(length(pjs)): $(pjs[p].name)")

    # Get construction time range for this project's scale
    construction_time_range = construction_time_ranges[pjs[p].scale]

    # Run Shapley sensitivity analysis with QUIET MODE enabled
    shapley_results = shapley_sensitivity_index(
        opt_scaling, n, wacc, electricity_price_mean, pjs[p];
        construction_time_range=construction_time_range,
        quiet=true  # QUIET MODE: Suppresses coalition-level logging
    )

    # Extract Shapley effects
    shapley_npv_results.res = [
        shapley_results.sh_npv.wacc,
        shapley_results.sh_npv.construction_time,
        shapley_results.sh_npv.loadfactor,
        shapley_results.sh_npv.investment
    ]
    rename!(shapley_npv_results, :res => pjs[p].name)

    shapley_lcoe_results.res = [
        shapley_results.sh_lcoe.wacc,
        shapley_results.sh_lcoe.construction_time,
        shapley_results.sh_lcoe.loadfactor,
        shapley_results.sh_lcoe.investment
    ]
    rename!(shapley_lcoe_results, :res => pjs[p].name)

    @info("Shapley sensitivity results",
          pj.name=pjs[p].name,
          pj.type=pjs[p].type,
          Sh_NPV=shapley_results.sh_npv,
          Sh_LCOE=shapley_results.sh_lcoe)

    @info("Completed $(p)/$(length(pjs)) reactors")

    # Force garbage collection between reactors to prevent memory buildup
    GC.gc()
end

# Output Shapley results
CSV.write("$outputpath/shapley-npv_results-$opt_scaling.csv", shapley_npv_results)
CSV.write("$outputpath/shapley-lcoe_results-$opt_scaling.csv", shapley_lcoe_results)

@info("="^80)
@info("SHAPLEY SENSITIVITY ANALYSIS COMPLETE")
@info("="^80)
@info("Results saved to:")
@info("  - shapley-npv_results-$opt_scaling.csv")
@info("  - shapley-lcoe_results-$opt_scaling.csv")
@info("")
@info("All analysis steps complete!")
