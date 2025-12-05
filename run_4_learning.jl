##### Run Learning Curve Analysis #####
# Part 4 of 4-part workflow
# Runtime: ~3-4 hours (30 N × 3 LR × 41 reactors × 100k samples)
# Outputs: learning-lcoe_curves-*.csv

using Pkg
Pkg.activate(pwd())
using Statistics
using CSV
using DataFrames

inputpath = "_input"
outputpath = "_output"

@info("="^80)
@info("STEP 4/4: LEARNING CURVE ANALYSIS")
@info("="^80)

##### load functions #####
@info("Loading functions")
include("functions.jl")

##### load project data #####
@info("Loading data")
include("data.jl")

##### simulation parameters #####

# number of Monte Carlo runs per scenario
n = Int64(100_000)  # 100k samples per learning scenario

# wholesale electricity price [EUR/MWh] - fixed at mean value
electricity_price_mean = mean([52.2, 95.8])

# weighted average cost of capital (WACC), lower and upper bound
wacc = [0.04, 0.10]

# construction time ranges by scale [years]
construction_time_ranges = Dict(
    "Micro" => [3, 8],   # Thesis: 3-8 years, triangular mode 5 technology → wider range (3-7 years)
    "SMR"   => [3, 8],   # Thesis: 3-8 years, triangular mode 5 technology → wider range (3-7 years)
    "Large" => [5, 13]   # Thesis: 5-13 years, triangular mode 8 data: Korea 5-6 yrs, US/Europe 7-12 yrs
)

# scaling parameter
scaling = [0.4, 0.7]  # Thesis: uniform distribution

##### Configuration: Load centralized settings #####
include("config.jl")



##### LEARNING CURVE PARAMETERS #####

# Cumulative units to simulate
N_units = [1, 2, 4, 6, 8, 12, 16, 20, 24, 30]

# Learning rates to test
learning_rates = [0.05, 0.10, 0.15]  # Pessimistic, base, optimistic
lr_labels = Dict(0.05 => "Pessimistic", 0.10 => "Base", 0.15 => "Optimistic")

# FOAK parameters
kappa = 1.0  # FOAK at baseline (no premium)
floor_m = nothing  # Unlimited learning (no floor)

@info("Learning curve configuration:")
@info("  N_units: $N_units")
@info("  Learning rates: $learning_rates")
@info("  FOAK premium (kappa): $kappa")
@info("  Learning floor: $floor_m")

##### RUN LEARNING CURVE ANALYSIS #####

@info("Starting learning curve analysis")
@info("Reactors: $(length(pjs))")
@info("Scenarios per reactor: $(length(N_units) * length(learning_rates))")
@info("Samples per scenario: $n")
@info("Total simulations: $(n * length(pjs) * length(N_units) * length(learning_rates))")

# Initialize results DataFrame
learning_results = DataFrame(
    reactor = String[],
    scale = String[],
    type = String[],
    N_unit = Int[],
    LR = Float64[],
    LR_label = String[],
    LCOE_mean = Float64[],
    LCOE_std = Float64[],
    LCOE_median = Float64[],
    LCOE_p10 = Float64[],
    LCOE_p90 = Float64[],
    LCOE_q25 = Float64[],
    LCOE_q75 = Float64[]
)

# Track progress
total_scenarios = length(pjs) * length(N_units) * length(learning_rates)
scenario_count = 0
start_time = time()

# Nested loops: reactor → N_unit → learning rate
for (pj_idx, pj) in enumerate(pjs)
    @info("="^60)
    @info("Reactor $(pj_idx)/$(length(pjs)): $(pj.name) ($(pj.scale) $(pj.type))")

    # Get construction time range for this reactor's scale
    construction_time_range = construction_time_ranges[pj.scale]

    for N in N_units
        for LR in learning_rates
            global scenario_count
            scenario_count += 1
            elapsed = time() - start_time
            avg_time_per_scenario = scenario_count > 1 ? elapsed / (scenario_count - 1) : 0
            remaining_scenarios = total_scenarios - scenario_count
            eta_seconds = remaining_scenarios * avg_time_per_scenario
            eta_minutes = eta_seconds / 60

            @info("  Scenario $scenario_count/$total_scenarios: N=$N, LR=$LR ($(lr_labels[LR])) | ETA: $(round(eta_minutes, digits=1)) min")

            # Generate random variables with learning
            rand_vars = gen_rand_vars(
                opt_scaling, n, wacc, electricity_price_mean, pj;
                apply_learning = true,
                N_unit = N,
                LR = LR,
                kappa = kappa,
                floor_m = floor_m,
                construction_time_range = construction_time_range
            )

            # Run investment simulation
            results = investment_simulation(pj, rand_vars)
            lcoe_values = vec(results[2])

            # Calculate statistics
            push!(learning_results, [
                pj.name,
                pj.scale,
                pj.type,
                N,
                LR,
                lr_labels[LR],
                mean(lcoe_values),
                std(lcoe_values),
                median(lcoe_values),
                quantile(lcoe_values, 0.10),
                quantile(lcoe_values, 0.90),
                quantile(lcoe_values, 0.25),
                quantile(lcoe_values, 0.75)
            ])
        end
    end
end

total_time = (time() - start_time) / 60
@info("="^80)
@info("Learning curve analysis complete!")
@info("Total runtime: $(round(total_time, digits=1)) minutes")
@info("Results: $(nrow(learning_results)) scenarios")

##### SAVE RESULTS #####

output_file = "$outputpath/learning-lcoe_curves-$opt_scaling.csv"
CSV.write(output_file, learning_results)
@info("Saved: $output_file")

@info("="^80)
@info("LEARNING CURVE ANALYSIS FINISHED")
@info("="^80)
