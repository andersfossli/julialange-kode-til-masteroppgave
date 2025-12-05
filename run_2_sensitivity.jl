##### Run ONLY Sobol Sensitivity Index Analysis #####
# Part 2 of 3-part workflow
# Runtime: ~30-60 minutes with 100k samples
# Outputs: si-npv_results-*.csv, si-lcoe_results-*.csv

using Pkg
Pkg.activate(pwd())
using Statistics

inputpath = "_input"
outputpath = "_output"

@info("="^80)
@info("STEP 2/3: SOBOL SENSITIVITY INDEX ANALYSIS")
@info("="^80)

##### load functions #####
@info("Loading functions")
include("functions.jl")

##### load project data #####
@info("Loading data")
include("data.jl")

##### simulation parameters #####

# number of Monte Carlo runs
n = Int64(100000)  # PRODUCTION: 100k samples for robust estimates

# wholesale electricity price [EUR/MWh]
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

##### RUN SENSITIVITY INDEX ANALYSIS ONLY #####

@info("Starting Sobol sensitivity index analysis")
@info("Reactors: $(length(pjs))")
@info("Samples per reactor: $n")

# initialize results variables
si_npv_results = DataFrame()
si_npv_results.si = ["S", "S", "S", "S", "ST", "ST", "ST", "ST"]
si_npv_results.var = ["wacc", "construction_time", "capacity factor", "scaling", "wacc", "construction_time", "capacity factor", "scaling"]
si_lcoe_results = DataFrame()
si_lcoe_results.si = ["S", "S", "S", "S", "ST", "ST", "ST", "ST"]
si_lcoe_results.var = ["wacc", "construction_time", "capacity factor", "scaling", "wacc", "construction_time", "capacity factor", "scaling"]

# run sensitivity analysis for all projects
for p in eachindex(pjs)
    @info("Running sensitivity analysis for reactor $(p)/$(length(pjs)): $(pjs[p].name)")

    # get construction time range for this project's scale
    construction_time_range = construction_time_ranges[pjs[p].scale]

    si_results = sensitivity_index(opt_scaling, n, wacc, electricity_price_mean, pjs[p];
                                   construction_time_range=construction_time_range)

    si_npv_results.res = vcat(collect(si_results[1]), collect(si_results[2]))
    rename!(si_npv_results, :res => pjs[p].name)
    si_lcoe_results.res = vcat(collect(si_results[3]), collect(si_results[4]))
    rename!(si_lcoe_results, :res => pjs[p].name)
end

# output
CSV.write("$outputpath/si-npv_results-$opt_scaling.csv", si_npv_results)
CSV.write("$outputpath/si-lcoe_results-$opt_scaling.csv", si_lcoe_results)

@info("="^80)
@info("SENSITIVITY INDEX ANALYSIS COMPLETE")
@info("="^80)
@info("Results saved to:")
@info("  - si-npv_results-$opt_scaling.csv")
@info("  - si-lcoe_results-$opt_scaling.csv")
@info("")
@info("Next step: Run run_3_shapley.jl for Shapley sensitivity analysis")
