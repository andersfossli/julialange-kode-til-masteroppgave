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

##### Define learning scenarios #####
# Each scenario tuple contains: (apply_learning, N_unit, LR, kappa, floor_m, tag)
# - apply_learning: true/false to enable learning curve
# - N_unit: number of units built (experience level)
# - LR: learning rate (e.g., 0.10 = 10% cost reduction per doubling)
# - kappa: FOAK premium multiplier (e.g., 1.20 = 20% above SOAK)
# - floor_m: minimum multiplier (e.g., 1.0 = SOAK baseline), use nothing for no floor
# - tag: label for output files

learning_cases = [
    # Baseline: no learning applied
    (false, 1, 0.0, 1.0, nothing, "baseline"),

    # Learning scenarios with LR=10%, FOAK premium κ=1.20 (20% above SOAK), floor at 1.0 (SOAK)
    (true,  1, 0.10, 1.20, 1.0, "LR10_N1_k120"),    # FOAK: 20% premium
    (true,  2, 0.10, 1.20, 1.0, "LR10_N2_k120"),    # 2nd unit: ~8% premium
    (true,  4, 0.10, 1.20, 1.0, "LR10_N4_k120"),    # 4th unit: hits SOAK floor
    (true,  6, 0.10, 1.20, 1.0, "LR10_N6_k120"),    # 6th unit: at SOAK

    # Optional: Higher learning rate scenarios (LR=15%)
    # (true,  1, 0.15, 1.20, 1.0, "LR15_N1_k120"),
    # (true,  2, 0.15, 1.20, 1.0, "LR15_N2_k120"),
    # (true,  4, 0.15, 1.20, 1.0, "LR15_N4_k120"),

    # Optional: No floor scenarios (learning continues below SOAK)
    # (true,  8, 0.10, 1.20, nothing, "LR10_N8_nofloor"),
]

@info("Running $(length(learning_cases)) learning scenarios")

##### Run simulation for each learning scenario #####

for (applyL, Nunit, LR, kappa, floor_m, tag) in learning_cases
    @info("=" ^ 80)
    @info("Learning scenario: $tag")
    @info("  apply_learning=$applyL, N_unit=$Nunit, LR=$LR, κ=$kappa, floor=$floor_m")
    @info("=" ^ 80)

    # Initialize result variables for this scenario
    npv_results = DataFrame()
    lcoe_results = DataFrame()

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
    npv_summary = describe(npv_results, :all)
    lcoe_summary = describe(lcoe_results, :all)

    # Save outputs with scenario tag
    CSV.write("$outputpath/mcs-npv_results-$opt_scaling-$tag.csv", npv_results)
    CSV.write("$outputpath/mcs-npv_summary-$opt_scaling-$tag.csv", npv_summary[!,1:8])
    CSV.write("$outputpath/mcs-lcoe_results-$opt_scaling-$tag.csv", lcoe_results)
    CSV.write("$outputpath/mcs-lcoe_summary-$opt_scaling-$tag.csv", lcoe_summary[!,1:8])

    @info("Saved results for scenario: $tag")
end

@info("=" ^ 80)
@info("All learning scenarios complete!")
@info("Output files saved to: $outputpath/")
@info("=" ^ 80)
