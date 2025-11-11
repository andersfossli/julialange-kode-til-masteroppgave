##### CENTRALIZED CONFIGURATION #####
# This file contains settings used by all simulation and plotting scripts
# Change parameters here to apply to all workflows

##### SCALING METHOD CONFIGURATION #####
# scaling options: 1=manufacturer, 2=roulstone, 3=rothwell, 4=uniform, 5=carelli
opts_scaling = ["manufacturer", "roulstone", "rothwell", "uniform", "carelli"]

# THESIS DEFAULT: Select scaling method (change here to apply everywhere)
local_scaling_index = 2  # 2 = roulstone (thesis default)

# Determine scaling option (supports both cluster and interactive mode)
if @isdefined(par_job) == true
    opt_scaling = opts_scaling[par_job]
    @info("Cluster job mode: using scaling option $opt_scaling (index $par_job)")
else
    opt_scaling = opts_scaling[local_scaling_index]
    @info("Interactive mode: using scaling option $opt_scaling (index $local_scaling_index)")
end

@info("="^80)
@info("CONFIGURATION LOADED: Scaling method = $opt_scaling")
@info("="^80)
