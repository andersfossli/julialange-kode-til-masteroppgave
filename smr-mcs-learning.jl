##### activate environment #####
using Pkg
Pkg.activate(pwd())
using Statistics
using CSV
using DataFrames

inputpath = "_input"
outputpath = "_output"

##### load functions #####
@info("Loading functions")
include("functions.jl");

##### load project data #####
@info("loading data")
include("data.jl");

##### further simulation data #####

    # number of Monte Carlo runs
    n = Int64(10000);

    # wholesale electricity price [USD/MWh], lower and upper bound of rand variable
    electricity_price = [52.2, 95.8];

    # weighted average cost of capital (WACC), lower and upper bound of rand variable
    wacc = [0.04, 0.10];

    # scaling
        # scaling options
        opts_scaling = ["manufacturer", "roulstone", "rothwell", "uniform"];
        # scaling parameter, lower and upper bound of random variable
        scaling = [0.20, 0.75];

    # choose scaling option
    if @isdefined(par_job) == true
        # read scaling option from job script parameter
        opt_scaling = opts_scaling[par_job];
        @info("using scaling option $opt_scaling")
    else
        # define scaling option locally
        opt_scaling = opts_scaling[2];  # default: roulstone
    end

##### Configuration #####
# Set to true to save NPV files in addition to LCOE files
save_npv_files = false  # Default: only save LCOE to reduce clutter

##### Define learning scenarios #####
# Each scenario tuple contains: (apply_learning, N_unit, LR, kappa, floor_m, tag)
# - apply_learning: true/false to enable learning curve
# - N_unit: number of units built (experience level)
# - LR: learning rate (e.g., 0.10 = 10% cost reduction per doubling)
# - kappa: FOAK premium multiplier (e.g., 1.20 = 20% above SOAK)
# - floor_m: minimum multiplier (e.g., 1.0 = SOAK baseline), use nothing for no floor
# - tag: label for output files

# Simplified default scenarios (3 cases: baseline, FOAK, SOAK)
# This produces 6 files total (3 scenarios × 2 files each)
learning_cases = [
    # Baseline: no learning applied (reference case)
    (false, 1, 0.0, 1.0, nothing, "baseline"),

    # FOAK: First unit with 20% premium
    (true,  1, 0.10, 1.20, 1.0, "FOAK"),

    # SOAK: Fourth unit (learned down to baseline)
    (true,  4, 0.10, 1.20, 1.0, "SOAK"),
]

# Extended scenarios - uncomment to run full learning curve analysis
# WARNING: This will create 24 files (6 scenarios × 4 files each)
# learning_cases = [
#     (false, 1, 0.0, 1.0, nothing, "baseline"),
#     (true,  1, 0.10, 1.20, 1.0, "LR10_N1_k120"),    # FOAK: 20% premium
#     (true,  2, 0.10, 1.20, 1.0, "LR10_N2_k120"),    # 2nd unit: ~8% premium
#     (true,  4, 0.10, 1.20, 1.0, "LR10_N4_k120"),    # 4th unit: hits SOAK floor
#     (true,  6, 0.10, 1.20, 1.0, "LR10_N6_k120"),    # 6th unit: at SOAK
#     (true,  8, 0.10, 1.20, 1.0, "LR10_N8_k120"),    # 8th unit: at SOAK
# ]

# Alternative: Higher learning rate scenarios (LR=15%)
# learning_cases = [
#     (false, 1, 0.0, 1.0, nothing, "baseline"),
#     (true,  1, 0.15, 1.20, 1.0, "LR15_FOAK"),
#     (true,  4, 0.15, 1.20, 1.0, "LR15_SOAK"),
# ]

@info("Running $(length(learning_cases)) learning scenarios")

##### Run simulation for each learning scenario #####

for (applyL, Nunit, LR, kappa, floor_m, tag) in learning_cases
    @info("=" ^ 80)
    @info("Learning scenario: $tag")
    @info("  apply_learning=$applyL, N_unit=$Nunit, LR=$LR, κ=$kappa, floor=$floor_m")
    @info("=" ^ 80)

    # Initialize result variables for this scenario
    local npv_results = DataFrame()
    local lcoe_results = DataFrame()

    # Run simulation for all projects
    for p in eachindex(pjs)
        @info("Running simulation for reactor $(pjs[p].name) ($(pjs[p].scale))")

        # Generate random variables with learning parameters
        rand_vars = gen_rand_vars(opt_scaling, n, wacc, electricity_price, pjs[p];
                                  apply_learning=applyL,
                                  N_unit=Nunit,
                                  LR=LR,
                                  kappa=kappa,
                                  floor_m=floor_m)

        # Run Monte Carlo simulation
        results = investment_simulation(pjs[p], rand_vars)

        # Normalize NPV to plant capacity [USD/MW]
        npv_results.res = vec(results[1] / pjs[p].plant_capacity)
        rename!(npv_results, :res => pjs[p].name)

        lcoe_results.res = vec(results[2])
        rename!(lcoe_results, :res => pjs[p].name)
    end

    # Calculate summary statistics
    local lcoe_summary = describe(lcoe_results, :all)

    # Save LCOE outputs (always)
    CSV.write("$outputpath/mcs-lcoe_results-$opt_scaling-$tag.csv", lcoe_results)
    CSV.write("$outputpath/mcs-lcoe_summary-$opt_scaling-$tag.csv", lcoe_summary[!,1:8])

    # Save NPV outputs (optional - controlled by save_npv_files flag)
    if save_npv_files
        local npv_summary = describe(npv_results, :all)
        CSV.write("$outputpath/mcs-npv_results-$opt_scaling-$tag.csv", npv_results)
        CSV.write("$outputpath/mcs-npv_summary-$opt_scaling-$tag.csv", npv_summary[!,1:8])
        @info("Saved LCOE and NPV results for scenario: $tag")
    else
        @info("Saved LCOE results for scenario: $tag (NPV files not saved, set save_npv_files=true to enable)")
    end
end

@info("=" ^ 80)
@info("All learning scenarios complete!")
num_files = length(learning_cases) * (save_npv_files ? 4 : 2)
@info("Created $num_files output files in: $outputpath/")
@info("  - LCOE results and summary for each scenario")
if save_npv_files
    @info("  - NPV results and summary for each scenario")
end
@info("=" ^ 80)
