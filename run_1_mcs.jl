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

# wholesale electricity price [EUR/MWh] - now fixed at mean value
electricity_price_mean = mean([52.2, 95.8])

# weighted average cost of capital (WACC), lower and upper bound
wacc = [0.04, 0.10]

# construction time ranges by scale [years]
construction_time_ranges = Dict(
    "Micro" => [3, 8],   # Thesis: 3-8 years, triangular mode 5
    "SMR"   => [3, 7],   # Thesis: 3-7 years, triangular mode 5
    "Large" => [5, 13]   # Thesis: 5-13 years, triangular mode 8
)

# scaling parameter (thesis: uniform distribution)
scaling = [0.4, 0.7]

##### Configuration: Load centralized settings #####
include("config.jl")

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

# initialize component breakdown dataframes
lcoe_occ_results = DataFrame()
lcoe_idc_results = DataFrame()
lcoe_fixed_om_results = DataFrame()
lcoe_variable_om_fuel_results = DataFrame()

# run simulation for all projects
for p in eachindex(pjs)
    @info("Running MCS for reactor $(p)/$(length(pjs)): $(pjs[p].name)")

    # get construction time range for this project's scale
    construction_time_range = construction_time_ranges[pjs[p].scale]

    # generate random variables
    rand_vars = gen_rand_vars(opt_scaling, n, wacc, electricity_price_mean, pjs[p];
                              construction_time_range=construction_time_range)

    # run Monte Carlo simulation with component tracking
    disc_res = mc_run(n, pjs[p], rand_vars)

    # calculate NPV and LCOE with component breakdown
    npv_lcoe_res = npv_lcoe(disc_res, decompose=true)

    # normalize NPV to plant capacity [EUR/MW]
    npv_results.res = vec(npv_lcoe_res.npv / pjs[p].plant_capacity)
    rename!(npv_results, :res => pjs[p].name)
    lcoe_results.res = vec(npv_lcoe_res.lcoe)
    rename!(lcoe_results, :res => pjs[p].name)

    # Save component breakdowns
    lcoe_occ_results.res = vec(npv_lcoe_res.lcoe_occ)
    rename!(lcoe_occ_results, :res => pjs[p].name)
    lcoe_idc_results.res = vec(npv_lcoe_res.lcoe_idc)
    rename!(lcoe_idc_results, :res => pjs[p].name)
    lcoe_fixed_om_results.res = vec(npv_lcoe_res.lcoe_fixed_om)
    rename!(lcoe_fixed_om_results, :res => pjs[p].name)
    lcoe_variable_om_fuel_results.res = vec(npv_lcoe_res.lcoe_variable_om_fuel)
    rename!(lcoe_variable_om_fuel_results, :res => pjs[p].name)

    # Save WACC and investment values 
    wacc_values.res = vec(rand_vars.wacc)
    rename!(wacc_values, :res => pjs[p].name)
    investment_values.res = vec(rand_vars.investment)
    rename!(investment_values, :res => pjs[p].name)
end

# summary statistics
npv_summary = describe(npv_results, :all)
lcoe_summary = describe(lcoe_results, :all)

# output - main results
CSV.write("$outputpath/mcs-npv_results-$opt_scaling.csv", npv_results)
CSV.write("$outputpath/mcs-npv_summary-$opt_scaling.csv", npv_summary[!,1:8])
CSV.write("$outputpath/mcs-lcoe_results-$opt_scaling.csv", lcoe_results)
CSV.write("$outputpath/mcs-lcoe_summary-$opt_scaling.csv", lcoe_summary[!,1:8])
CSV.write("$outputpath/mcs-wacc_values-$opt_scaling.csv", wacc_values)
CSV.write("$outputpath/mcs-investment_values-$opt_scaling.csv", investment_values)

# output - LCOE component breakdown
CSV.write("$outputpath/mcs-lcoe_occ-$opt_scaling.csv", lcoe_occ_results)
CSV.write("$outputpath/mcs-lcoe_idc-$opt_scaling.csv", lcoe_idc_results)
CSV.write("$outputpath/mcs-lcoe_fixed_om-$opt_scaling.csv", lcoe_fixed_om_results)
CSV.write("$outputpath/mcs-lcoe_variable_om_fuel-$opt_scaling.csv", lcoe_variable_om_fuel_results)

@info("="^80)
@info("MONTE CARLO SIMULATION COMPLETE")
@info("="^80)
@info("Results saved to:")
@info("  - mcs-npv_results-$opt_scaling.csv")
@info("  - mcs-lcoe_results-$opt_scaling.csv")
@info("  - mcs-*_summary-$opt_scaling.csv")
@info("  - mcs-wacc_values-$opt_scaling.csv")
@info("  - mcs-investment_values-$opt_scaling.csv")
@info("  - mcs-lcoe_occ-$opt_scaling.csv (NEW: OCC component)")
@info("  - mcs-lcoe_idc-$opt_scaling.csv (NEW: IDC component)")
@info("  - mcs-lcoe_fixed_om-$opt_scaling.csv (NEW: Fixed O&M component)")
@info("  - mcs-lcoe_variable_om_fuel-$opt_scaling.csv (NEW: Variable O&M+Fuel component)")
@info("")
@info("Next step: Run run_2_sensitivity.jl for Sobol sensitivity indices")
