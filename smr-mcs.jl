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
        opt_scaling = opts_scaling[2];
    end

##### run simulation #####

@info("running simulation")
include("run_simulation.jl")