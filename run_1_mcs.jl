##### Run ONLY Monte Carlo Simulation #####
# Part 1 of 3-part workflow
# Runtime: ~1-2 hours with 100k samples
# Outputs: mcs-*_results-*.csv, mcs-*_summary-*.csv

using Pkg
Pkg.activate(pwd())
using Statistics

inputpath = "_input"
outputpath = "_output"

@info("="^80)
@info("STEP 1/3: MONTE CARLO SIMULATION")
@info("="^80)

##### load functions #####
@info("Loading functions")
include("functions.jl")

##### load project data #####
@info("Loading data")
include("data.jl")

##### simulation parameters #####

# number of Monte Carlo runs
n = Int64(1_000_000)  # PRODUCTION: 1M samples for publication-quality estimates (matches Uncertainties report)

# wholesale electricity price [USD/MWh] - now fixed at mean value
electricity_price_mean = mean([52.2, 95.8])

# weighted average cost of capital (WACC), lower and upper bound
wacc = [0.04, 0.10]

# construction time ranges by scale [years]
construction_time_ranges = Dict(
    "Micro" => [3, 7],   # Unproven technology → wider range (3-7 years)
    "SMR"   => [3, 7],   # Unproven technology → wider range (3-7 years)
    "Large" => [5, 12]   # Historical data: Korea 5-6 yrs, US/Europe 7-12 yrs
)

# scaling parameter
scaling = [0.20, 0.75]

# scaling options: 1=manufacturer, 2=roulstone, 3=rothwell, 4=uniform, 5=carelli
opts_scaling = ["manufacturer", "roulstone", "rothwell", "uniform", "carelli"]

# CONFIGURATION: Select scaling method
local_scaling_index = 3  # 3 = rothwell

# choose scaling option
if @isdefined(par_job) == true
    opt_scaling = opts_scaling[par_job]
    @info("Cluster job mode: using scaling option $opt_scaling (index $par_job)")
else
    opt_scaling = opts_scaling[local_scaling_index]
    @info("Interactive mode: using scaling option $opt_scaling (index $local_scaling_index)")
end

##### RUN MONTE CARLO SIMULATION ONLY #####

@info("Starting Monte Carlo simulation")
@info("Reactors: $(length(pjs))")
@info("Samples per reactor: $n")
@info("Total simulations: $(n * length(pjs))")

# initialize result variables
npv_results = DataFrame()
lcoe_results = DataFrame()
wacc_values = DataFrame()
investment_values = DataFrame()

# run simulation for all projects
for p in eachindex(pjs)
    @info("Running MCS for reactor $(p)/$(length(pjs)): $(pjs[p].name)")

    # get construction time range for this project's scale
    construction_time_range = construction_time_ranges[pjs[p].scale]

    # generate random variables
    rand_vars = gen_rand_vars(opt_scaling, n, wacc, electricity_price_mean, pjs[p];
                              construction_time_range=construction_time_range)

    # run Monte Carlo simulation
    results = investment_simulation(pjs[p], rand_vars)

    # normalize NPV to plant capacity [USD/MW]
    npv_results.res = vec(results[1] / pjs[p].plant_capacity)
    rename!(npv_results, :res => pjs[p].name)
    lcoe_results.res = vec(results[2])
    rename!(lcoe_results, :res => pjs[p].name)

    # Save WACC and investment values
    wacc_values.res = vec(rand_vars.wacc)
    rename!(wacc_values, :res => pjs[p].name)
    investment_values.res = vec(rand_vars.investment)
    rename!(investment_values, :res => pjs[p].name)
end

# summary statistics
npv_summary = describe(npv_results, :all)
lcoe_summary = describe(lcoe_results, :all)

# output
CSV.write("$outputpath/mcs-npv_results-$opt_scaling.csv", npv_results)
CSV.write("$outputpath/mcs-npv_summary-$opt_scaling.csv", npv_summary[!,1:8])
CSV.write("$outputpath/mcs-lcoe_results-$opt_scaling.csv", lcoe_results)
CSV.write("$outputpath/mcs-lcoe_summary-$opt_scaling.csv", lcoe_summary[!,1:8])
CSV.write("$outputpath/mcs-wacc_values-$opt_scaling.csv", wacc_values)
CSV.write("$outputpath/mcs-investment_values-$opt_scaling.csv", investment_values)

@info("="^80)
@info("MONTE CARLO SIMULATION COMPLETE")
@info("="^80)
@info("Results saved to:")
@info("  - mcs-npv_results-$opt_scaling.csv")
@info("  - mcs-lcoe_results-$opt_scaling.csv")
@info("  - mcs-*_summary-$opt_scaling.csv")
@info("  - mcs-wacc_values-$opt_scaling.csv")
@info("  - mcs-investment_values-$opt_scaling.csv")
@info("")
@info("Next step: Run run_2_sensitivity.jl for Sobol sensitivity indices")
