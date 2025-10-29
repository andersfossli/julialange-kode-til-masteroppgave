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

# Default scenarios: baseline + three learning rate scenarios
# κ = 1.0 (FOAK-anchored, no premium), floor = nothing (unlimited learning)
# This produces 32 files (16 scenarios × 2 LCOE files each)
learning_cases = [
    # Baseline: no learning applied (reference SOAK case)
    (false, 1, 0.00, 1.00, nothing,  "baseline"),

    # Conservative: LR=5% (slow learning)
    (true,  1, 0.05, 1.00, nothing,   "LR05_N1_k100"),    # FOAK: m=1.00
    (true,  2, 0.05, 1.00, nothing,   "LR05_N2_k100"),    # 2nd: m=0.95
    (true,  4, 0.05, 1.00, nothing,   "LR05_N4_k100"),    # 4th: m=0.90
    (true,  6, 0.05, 1.00, nothing,   "LR05_N6_k100"),    # 6th: m=0.87
    (true, 12, 0.05, 1.00, nothing,   "LR05_N12_k100"),   # 12th: m=0.83

    # Base: LR=10% (moderate learning)
    (true,  1, 0.10, 1.00, nothing,   "LR10_N1_k100"),    # FOAK: m=1.00
    (true,  2, 0.10, 1.00, nothing,   "LR10_N2_k100"),    # 2nd: m=0.90
    (true,  4, 0.10, 1.00, nothing,   "LR10_N4_k100"),    # 4th: m=0.81
    (true,  6, 0.10, 1.00, nothing,   "LR10_N6_k100"),    # 6th: m=0.76
    (true, 12, 0.10, 1.00, nothing,   "LR10_N12_k100"),   # 12th: m=0.68

    # Optimistic: LR=15% (fast learning)
    (true,  1, 0.15, 1.00, nothing,   "LR15_N1_k100"),    # FOAK: m=1.00
    (true,  2, 0.15, 1.00, nothing,   "LR15_N2_k100"),    # 2nd: m=0.85
    (true,  4, 0.15, 1.00, nothing,   "LR15_N4_k100"),    # 4th: m=0.72
    (true,  6, 0.15, 1.00, nothing,   "LR15_N6_k100"),    # 6th: m=0.66
    (true, 12, 0.15, 1.00, nothing,   "LR15_N12_k100"),   # 12th: m=0.56
]

# Alternative: Higher learning rate scenarios (LR=15%)
# Uncomment to explore more aggressive learning
# learning_cases = [
#     (false, 1, 0.00, 1.00, nothing, "baseline"),
#     (true,  1, 0.15, 1.20, 1.00, "LR15_N1_k120"),
#     (true,  2, 0.15, 1.20, 1.00, "LR15_N2_k120"),
#     (true,  4, 0.15, 1.20, 1.00, "LR15_N4_k120"),
#     (true,  6, 0.15, 1.20, 1.00, "LR15_N6_k120"),
#     (true, 12, 0.15, 1.20, 1.00, "LR15_N12_k120"),
# ]

# Simplified scenarios (if you want fewer files - 6 total)
# learning_cases = [
#     (false, 1, 0.0, 1.0, nothing, "baseline"),
#     (true,  1, 0.10, 1.20, 1.0, "FOAK"),
#     (true,  4, 0.10, 1.20, 1.0, "SOAK"),
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
