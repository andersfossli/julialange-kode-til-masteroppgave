##### activate environment #####
using Pkg
Pkg.activate(pwd())
using Statistics

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

    # wholesale electricity price [USD/MWh] - now fixed at mean value
    # (electricity price doesn't affect LCOE calculation, so uncertainty removed)
    electricity_price_mean = mean([52.2, 95.8]);

    # weighted average cost of capital (WACC), lower and upper bound of rand variable
    wacc = [0.04, 0.10];

    # construction time ranges by scale [years]
    # Construction time now replaces electricity price as an uncertain parameter
    # (construction time affects LCOE through interest during construction)
    construction_time_ranges = Dict(
        "Micro" => [3, 7],   # Unproven technology → wider range (3-7 years)
        "SMR"   => [3, 7],   # Unproven technology → wider range (3-7 years)
        "Large" => [5, 12]   # Historical data: Korea 5-6 yrs, US/Europe 7-12 yrs
    );

    # scaling
        # scaling options: 1=manufacturer, 2=roulstone, 3=rothwell, 4=uniform, 5=carelli
        opts_scaling = ["manufacturer", "roulstone", "rothwell", "uniform", "carelli"];
        # scaling parameter, lower and upper bound of random variable
        # Note: carelli scaling doesn't use this parameter (uses fixed β=0.20, P_ref=1200 MWe)
        scaling = [0.20, 0.75];

    # CONFIGURATION: Select scaling method for interactive runs
    # For cluster jobs, this is overridden by: julia job.jl <index>
    local_scaling_index = 3;  # 1=manufacturer, 2=roulstone, 3=rothwell, 4=uniform, 5=carelli

    # choose scaling option
    if @isdefined(par_job) == true
        # Cluster job mode: read scaling option from job script parameter
        opt_scaling = opts_scaling[par_job];
        @info("Cluster job mode: using scaling option $opt_scaling (index $par_job)")
    else
        # Interactive mode: use local configuration
        opt_scaling = opts_scaling[local_scaling_index];
        @info("Interactive mode: using scaling option $opt_scaling (index $local_scaling_index)")
    end

##### run simulation #####

# NOTE: Learning is DISABLED in base Monte Carlo simulation
# - apply_learning=false (default)
# - apply_soak_discount=false (default)
# For learning scenarios with explicit LR/N parameters, use smr-mcs-learning.jl instead

@info("running simulation")
include("run_simulation.jl")