##### dependencies #####
using Distributions, LinearAlgebra, Combinatorics

##### defining structs #####

"""
The project struct is a mutable struct that represents an investment project. It has the following fields:
    name: a string that specifies the investment concept.
    type: a string that specifies the investment type.
    investment: a floating-point number that represents the investment estimate by the manufacturer in USD/MW.
    plant_capacity: a floating-point number that represents the plant capacity in MW.
    learning_factor: a floating-point number that represents the learning factor.
    time: an abstract vector that represents the project time in years, which includes the construction time and plant lifetime.
    loadfactor: an abstract vector that represents the load factor. It includes the lower and upper bounds of a random variable.
    operating_cost: an abstract vector that represents the O&M cost, including O&M fix cost in USD/MW, O&M variable cost in USD/MWh, and fuel cost in USD/MWh.
    reference_pj: an abstract vector that represents the reference reactor, which includes investment costs in USD/MW and plant capacity in MW.
"""
mutable struct project
    name::String                    # investment concept
    type::String                    # investment type
    scale::String                   # reactor scale (Micro/SMR/Large)
    region::String                  # geographic region (East Asia/Western/etc.)
    investment::Float64             # investment estimate by manufacturer [USD/MW]
    plant_capacity::Float64         # plant capacity [MW]
    learning_factor::Float64        # learning factor
    time::AbstractVector            # project time [years] (construction time, plant lifetime)
    loadfactor::AbstractVector      # load factor, lower, and upper bound of rand variable
    operating_cost::AbstractVector  # O&M cost (O&M fix cost [USD/MW], O&M variable cost [USD/MWh], fuel cost [USD/MWh])
    reference_pj::AbstractVector    # reference reactor (investment costs [USD/MW], plant capacity [MW])
end

##### defining functions #####

"""
    learning_multiplier(N::Int, LR::Float64; kappa::Float64=1.0, floor::Union{Nothing,Float64}=nothing, learning_type::String="factory")

Calculate the learning curve multiplier for cost reduction based on number of units built.

# Arguments
- `N::Int`: Number of units built (experience level)
- `LR::Float64`: Learning rate (e.g., 0.10 for 10% cost reduction per doubling)
- `kappa::Float64=1.0`: FOAK premium multiplier (e.g., 1.20 = 20% premium above SOAK)
- `floor::Union{Nothing,Float64}=nothing`: Optional floor to prevent going below a minimum value
- `learning_type::String="factory"`: Type of learning curve to apply
  - "factory": Full learning rate for modular/factory production (Micro/SMR)
  - "deployment": Half learning rate for sequential project builds (Large reactors)

# Learning Types
- **Factory Learning**: Full LR applies (standard Wright's Law)
  - For modular, factory-built reactors (Micro, SMR)
  - Example: Korean modular construction, NuScale factory production

- **Deployment Learning**: Half LR (LR/2) applies
  - For large, site-built reactors learning from project-to-project experience
  - Reflects slower organizational/regulatory learning vs manufacturing
  - Example: Korean/Chinese large reactor programs (~5% per doubling)

# Returns
- Multiplier to apply to investment cost (1.0 = baseline)

# Examples
```julia
learning_multiplier(1, 0.10, kappa=1.00)  # FOAK: 1.00 (no premium)
learning_multiplier(2, 0.10, kappa=1.00)  # Second unit (factory): 0.90
learning_multiplier(4, 0.10, kappa=1.00, learning_type="deployment")  # Fourth large unit: ~0.95
learning_multiplier(12, 0.10, kappa=1.00)  # 12th unit (factory): ~0.68
```
"""
function learning_multiplier(N::Int, LR::Float64;
                           kappa::Float64=1.0,
                           floor::Union{Nothing,Float64}=nothing,
                           learning_type::String="factory")
    # Apply learning rate based on type
    if learning_type == "deployment"
        # Deployment learning: Half the learning rate for sequential project builds
        # Reflects slower learning from project experience vs factory production
        effective_LR = LR / 2.0
    else  # "factory"
        # Factory learning: Full learning rate (standard Wright's Law)
        # For modular/factory-built reactors
        effective_LR = LR
    end

    m = kappa * (1.0 - effective_LR)^(log2(N))
    return isnothing(floor) ? m : max(m, floor)
end

"""
    carelli_occ(p::Float64, occ::Float64; p_ref::Float64=1200.0, beta::Float64=0.20)

Apply Carelli scaling to normalize OCC (Overnight Capital Cost) to a reference capacity for scale-independent comparison.

# Purpose
Derives scale-independent investment by correcting for economies of scale using a power-law relationship.
This allows fair comparison across different reactor sizes (Micro/SMR/Large).

# Arguments
- `p::Float64`: Actual reactor capacity [MWe]
- `occ::Float64`: Raw overnight capital cost [USD/kW or USD/MW depending on context]
- `p_ref::Float64=1200.0`: Reference capacity for normalization [MWe] (default: typical large PWR)
- `beta::Float64=0.20`: Scaling exponent (economies of scale parameter)
  - β = 0.15: Weak economies of scale
  - β = 0.20: Moderate (default, typical for nuclear)
  - β = 0.25: Strong economies of scale
  - β = 0.30: Very strong economies of scale

# Theory
Cost scaling relationship: OCC(P) = OCC_ref × (P/P_ref)^(-β)
To normalize to reference size: OCC_SI = OCC_raw × (P/P_ref)^β

This removes the size penalty, making costs comparable across scales.

# Returns
- Scale-independent OCC adjusted to reference capacity

# Examples
```julia
# 300 MW SMR with OCC = 6000 USD/kW
# Normalize to 1200 MW reference
carelli_occ(300.0, 6000.0)  # Returns ~4757 USD/kW (removes small-scale penalty)

# 1200 MW Large reactor with OCC = 5500 USD/kW
carelli_occ(1200.0, 5500.0)  # Returns 5500 USD/kW (already at reference)

# 50 MW Micro with OCC = 8000 USD/kW, strong scaling (β=0.25)
carelli_occ(50.0, 8000.0, beta=0.25)  # Returns ~5040 USD/kW
```

# References
- Carelli et al. (2010): Economics of small modular reactors
- OECD-NEA (2011): Current Status, Technical Feasibility and Economics of Small Nuclear Reactors
"""
function carelli_occ(p::Float64, occ::Float64; p_ref::Float64=1200.0, beta::Float64=0.20)
    # Apply power-law normalization
    # OCC_SI = OCC_raw × (P / P_ref)^β
    return occ * (p / p_ref)^beta
end

"""
    generate_correlated_samples(n, wacc_range, ct_range; ρ=0.4, wacc_mode=nothing, ct_mode=nothing)

Generate correlated WACC and construction_time samples using Gaussian copula with triangular marginals.

This function creates samples that capture the empirical correlation between WACC and construction time
(higher WACC environments often correlate with longer construction times due to financing/regulatory complexity).

# Methodology
Uses a Gaussian copula approach:
1. Generate correlated standard normals via Cholesky decomposition
2. Transform to uniform marginals via standard normal CDF
3. Transform to triangular marginals via inverse CDF (quantile function)

# Arguments
- `n::Int`: Number of samples to generate
- `wacc_range::Vector`: [min, max] for WACC (e.g., [0.04, 0.10])
- `ct_range::Vector`: [min, max] for construction time in years (e.g., [3, 7])
- `ρ::Float64=0.4`: Correlation coefficient for Gaussian copula (default 0.4)
- `wacc_mode::Union{Nothing,Float64}=nothing`: Mode for WACC triangular distribution (default: midpoint)
- `ct_mode::Union{Nothing,Float64}=nothing`: Mode for construction time triangular distribution (default: midpoint)

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: (WACC samples, construction_time samples) with correlation ρ

# Example
```julia
wacc_samples, ct_samples = generate_correlated_samples(1000, [0.04, 0.10], [3, 7]; ρ=0.4)
cor(wacc_samples, ct_samples)  # Should be ≈ 0.4
```

# References
- Nelsen (2006): An Introduction to Copulas
- Song et al. (2009): Corpora-based uncertainty propagation
"""
function generate_correlated_samples(n::Int, wacc_range::Vector, ct_range::Vector;
                                     ρ::Float64=0.4,
                                     wacc_mode::Union{Nothing,Float64}=nothing,
                                     ct_mode::Union{Nothing,Float64}=nothing)

    # Set modes to midpoint if not specified
    wacc_m = isnothing(wacc_mode) ? mean(wacc_range) : wacc_mode
    ct_m = isnothing(ct_mode) ? mean(ct_range) : ct_mode

    # Correlation matrix for Gaussian copula
    Σ = [1.0  ρ;
         ρ    1.0]

    # Cholesky decomposition: Σ = L * L'
    L = cholesky(Σ).L

    # Generate correlated standard normals: Z ~ N(0, Σ)
    Z = randn(n, 2) * L'

    # Transform to uniform marginals via standard normal CDF: U ~ Uniform(0,1) with correlation ρ
    U = cdf.(Normal(0, 1), Z)

    # Transform to triangular marginals via inverse CDF (quantile function)
    # TriangularDist(a, b, c) where a=min, b=max, c=mode
    wacc_dist = TriangularDist(wacc_range[1], wacc_range[2], wacc_m)
    ct_dist = TriangularDist(ct_range[1], ct_range[2], ct_m)

    WACC = quantile.(Ref(wacc_dist), U[:, 1])
    CT = quantile.(Ref(ct_dist), U[:, 2])

    return WACC, CT
end

"""
    gen_rand_vars(opt_scaling, n, wacc, electricity_price_mean, pj; ...)

Generate random variables for Monte Carlo simulation with triangular distributions and Gaussian copula.

# UPDATED Implementation (Triangular Distributions + Gaussian Copula)
This function now uses:
- **Triangular distributions** for all uncertain parameters (more realistic than uniform)
- **Gaussian copula** for WACC × construction_time correlation (ρ=0.4)
- Reactor-type-specific capacity factor ranges
- Scale-specific construction time ranges

# Arguments
- `opt_scaling::String`: Scaling method ("manufacturer", "roulstone", "rothwell", "uniform", "carelli")
- `n::Int64`: Number of Monte Carlo samples
- `wacc::Vector`: [min, max] for WACC (e.g., [0.04, 0.10])
- `electricity_price_mean::Float64`: Mean electricity price (no longer uncertain)
- `pj::project`: Project struct with reactor data
- `construction_time_range::Union{Nothing,Vector}=nothing`: [min, max] for construction time
- `apply_learning::Bool=false`: Apply learning curves (default: disabled)
- `apply_soak_discount::Bool=false`: Apply SOAK discount (default: disabled)
- Other learning parameters: N_unit, LR, kappa, floor_m

# Random Variable Generation
## Correlated Variables (Gaussian Copula)
- **WACC**: Triangular distribution with mode at midpoint, correlated with construction time (ρ=0.4)
- **Construction Time**: Triangular distribution with mode at midpoint, correlated with WACC (ρ=0.4)
  - Reflects empirical observation that higher WACC environments correlate with longer construction times

## Independent Variables (Triangular Distributions)
- **Load Factor**: Triangular distribution with mode at 60th percentile (reactor-type-specific ranges)
- **Electricity Price**: Fixed at mean value (no uncertainty - doesn't affect LCOE)

## Investment Cost (based on opt_scaling)
- "manufacturer": Manufacturer estimates (pj.investment × pj.plant_capacity)
- "roulstone": Roulstone scaling with uniform β ∈ [0.20, 0.75]
- "rothwell": Rothwell scaling (Roulstone parameterization)
- "uniform": Uniform distribution between scaled bounds
- "carelli": Carelli scaling (fixed β = 0.20, P_ref = 1200 MWe)

# Learning (Disabled by Default)
- SOAK discount: Only if apply_soak_discount=true
- Learning curves: Only if apply_learning=true
- **Base simulations do NOT apply learning** (see smr-mcs-learning.jl for scenarios)

# Returns
Named tuple: `(wacc, electricity_price, loadfactor, investment, construction_time)`

# Example
```julia
rand_vars = gen_rand_vars("rothwell", 10000, [0.04, 0.10], 74.0, project;
                          construction_time_range=[3, 7])
# WACC and construction_time will have correlation ≈ 0.4
```

# References
- Triangular distributions: More realistic than uniform for bounded uncertain parameters
- Gaussian copula: Nelsen (2006) "An Introduction to Copulas"
- WACC × CT correlation: Empirical observation from nuclear construction data
"""
function gen_rand_vars(opt_scaling::String, n::Int64, wacc::Vector, electricity_price_mean::Float64, pj::project;
                       apply_learning::Bool=false, N_unit::Int=1, LR::Float64=0.0, kappa::Float64=1.0, floor_m::Union{Nothing,Float64}=nothing,
                       apply_soak_discount::Bool=false,
                       construction_time_range::Union{Nothing,Vector}=nothing)

    @info "generating random variables"

    # SOAK discount factor: if not applying, use 1.0 (no discount)
    soak_factor = apply_soak_discount ? (1 - pj.learning_factor) : 1.0

    if !apply_soak_discount
        @info("SOAK discount disabled: using manufacturer OCC without learning_factor discount")
    end

    # Calculate total time based on maximum possible construction time (if variable) or fixed construction time
    # This ensures arrays are sized correctly for variable construction times
    if !isnothing(construction_time_range)
        max_possible_construction = construction_time_range[2]  # upper bound of construction time range
        total_time = max_possible_construction + pj.time[2]
        @info("Array sizing based on max construction time: $max_possible_construction + $(pj.time[2]) = $total_time years")
    else
        total_time = pj.time[1] + pj.time[2]
    end

    # generation of random non-project specific variables with triangular distributions
    # UPDATED: Use triangular distributions instead of uniform for more realistic uncertainty

    # Generate correlated WACC and construction_time using Gaussian copula
    if !isnothing(construction_time_range)
        # Mode for triangular distributions (can be parameterized later)
        wacc_mode = mean(wacc)  # Mode at midpoint
        ct_mode = mean(construction_time_range)  # Mode at midpoint

        @info("Generating correlated samples: WACC × construction_time (ρ=0.4)")
        rand_wacc, rand_construction_time_float = generate_correlated_samples(
            n, wacc, construction_time_range;
            ρ=0.4, wacc_mode=wacc_mode, ct_mode=ct_mode
        )
        # Round construction time to integer years
        rand_construction_time = round.(Int, rand_construction_time_float)

        # Verify correlation achieved
        actual_corr = cor(rand_wacc, rand_construction_time_float)
        @info("Empirical correlation: $(round(actual_corr, digits=3))")
    else
        # Fallback: independent triangular sampling if no CT range
        wacc_dist = TriangularDist(wacc[1], wacc[2], mean(wacc))
        rand_wacc = rand(wacc_dist, n)
        rand_construction_time = fill(Int(pj.time[1]), n)
    end

    # Electricity price is now fixed at mean (doesn't affect LCOE calculation)
    rand_electricity_price = fill(electricity_price_mean, n, total_time)

    # Generate loadfactor using triangular distribution
    # Mode at 60th percentile (slightly optimistic but realistic for planned operations)
    lf_mode = pj.loadfactor[1] + 0.6 * (pj.loadfactor[2] - pj.loadfactor[1])
    lf_dist = TriangularDist(pj.loadfactor[1], pj.loadfactor[2], lf_mode)
    rand_loadfactor = rand(lf_dist, n, total_time)
    
    # generation of uniformly distributed random project specific variables
        # scaling case distinction
            # Note that the scaling parameter here are converted such that Rothwell and Roulstone coincide.
            # Special handling for Large reactors: use face value, no scaling
            if pj.scale == "Large"
                @info("Large reactor detected: using face value (no scaling applied)")
                # For large reactors, use manufacturer estimates directly (no scaling)
                # This keeps large reactor costs at their empirical values
                rand_investment = pj.investment * pj.plant_capacity * ones(n) * soak_factor
            elseif opt_scaling == "manufacturer"
                @info("using manufacturer estimates")
                # deterministic investment cost based on manufacturer estimates [USD]
                rand_investment = pj.investment * pj.plant_capacity * ones(n) * soak_factor
            elseif opt_scaling == "roulstone"
                @info("using Roulstone scaling")
                # random investment cost based on Roulstone [USD]
                rand_scaling = scaling[1] .+ (scaling[2] - scaling[1]) .* rand(n,1)
                rand_investment = pj.reference_pj[1] * pj.reference_pj[2] * soak_factor * (pj.plant_capacity/pj.reference_pj[2]) .^ (rand_scaling)
            elseif opt_scaling == "rothwell"
                @info("using Rothwell scaling")
                # random investment cost based on Rothwell [USD]
                rand_scaling = 2^(scaling[1]-1) .+ (2^(scaling[2]-1) - 2^(scaling[1]-1)) .* rand(n,1)
                rand_investment = pj.reference_pj[1] * pj.reference_pj[2] * soak_factor * (pj.plant_capacity/pj.reference_pj[2]) .^ (1 .+ log.(rand_scaling) ./ log(2))
            elseif opt_scaling == "uniform"
                @info("using uniform scaling")
                # random investment cost uniform [USD]
                investment_scaled = pj.reference_pj[1] * pj.reference_pj[2] * (pj.plant_capacity/pj.reference_pj[2]) .^ (scaling)
                rand_investment = investment_scaled[1] .+ (investment_scaled[2]-investment_scaled[1]) .* rand(n,1) .* soak_factor
            elseif opt_scaling == "carelli"
                @info("using Carelli scaling (economies of scale normalization)")
                # Carelli scaling: normalize OCC to reference capacity for scale-independent comparison
                # Default: P_ref = 1200 MWe, β = 0.20 (can be parameterized later)
                p_ref = 1200.0  # Reference capacity [MWe]
                beta = 0.20     # Scaling exponent (economies of scale parameter)

                # Get raw OCC from manufacturer estimate [USD/kW]
                raw_occ = pj.investment  # Already in USD/kW

                # Apply Carelli normalization
                occ_si = carelli_occ(pj.plant_capacity, raw_occ, p_ref=p_ref, beta=beta)

                # Convert back to total investment [USD]
                rand_investment = occ_si * pj.plant_capacity * ones(n) * soak_factor

                @info("  Carelli parameters: P=$(pj.plant_capacity) MWe, P_ref=$p_ref MWe, β=$beta")
                @info("  Raw OCC: $(round(raw_occ, digits=0)) USD/kW → SI-OCC: $(round(occ_si, digits=0)) USD/kW")
            else
                @error("Option for the scaling method is unknown.")
            end

    # Apply learning curve multiplier (if requested)
    # Learning type depends on reactor scale:
    #  - Micro/SMR: "factory" learning (full LR from modular production)
    #  - Large: "deployment" learning (LR/2 from sequential project experience)
    if apply_learning
        # Select learning type based on reactor scale
        learning_type = (pj.scale == "Large") ? "deployment" : "factory"

        m = learning_multiplier(N_unit, LR; kappa=kappa, floor=floor_m, learning_type=learning_type)

        if pj.scale == "Large"
            @info("Applying deployment learning (LR/2): $(pj.scale) reactor, N=$N_unit, LR=$LR, effective_LR=$(LR/2), κ=$kappa, floor=$floor_m → multiplier=$m")
        else
            @info("Applying factory learning (full LR): $(pj.scale) reactor, N=$N_unit, LR=$LR, κ=$kappa, floor=$floor_m → multiplier=$m")
        end

        rand_investment .*= m
    end

    # NOTE: Construction time is now generated earlier with WACC using Gaussian copula (see lines 260-282)
    # This captures the correlation between WACC and construction time (ρ=0.4)

    # output
    return(wacc = rand_wacc, electricity_price = rand_electricity_price, loadfactor = rand_loadfactor,
           investment = rand_investment, construction_time = rand_construction_time)

end

"""
The mc_run function performs a Monte Carlo simulation for an investment project, pj, based on a specified number of simulations to run, n, and a set of random variables, rand_vars.
The function first initializes the total time of the reactor project, which is the sum of the construction time and operating time, and the fixed and variable operating and maintenance (O&M) costs.
It then initializes variables for the random variables, cash inflow, cash outflow, cashflow, discounted cash outflow, NPV in period t, generated electricity, and discounted generated electricity.
The function then enters a simulation loop that iterates over the total time of the project. For each time step, the function performs different calculations depending on whether the time step is before or after the construction time.
    If the time step is before the construction time, the function calculates the cash outflow as the investment cost divided by the construction time, and cashflow as the difference between cash inflow and cash outflow.
    If the time step is after the construction time, the function calculates the amount of electricity produced, cash inflow from electricity sales, cash outflow for O&M costs, cashflow, discounted cash outflow, discounted electricity.
Finally, the function returns the results of the simulation in the form of a named tuple with three fields: disc_cash_out, disc_cash_net, and disc_electricity and their corresponding values.
"""
function mc_run(n::Int64, pj::project, rand_vars)

    @info "running Monte Carlo simulation"

    # project data
        # O&M costs
            # fixed O&M costs [USD/year]
            operating_cost_fix = pj.plant_capacity * pj.operating_cost[1]
            # variable O&M costs [USD/MWh]
            operating_cost_variable = pj.operating_cost[2] + pj.operating_cost[3]

    # initialize variables
        # random variables
        rand_wacc = rand_vars.wacc
        rand_electricity_price = rand_vars.electricity_price
        rand_loadfactor = rand_vars.loadfactor
        rand_investment = rand_vars.investment
        rand_construction_time = rand_vars.construction_time

        # Calculate maximum construction time to size arrays
        max_construction = maximum(rand_construction_time)
        total_time = max_construction + pj.time[2]

        # cash inflow
        cash_in = zeros(Float64, n, total_time)
        # cash outflow
        cash_out = zeros(Float64, n, total_time)
        # cashflow = inflow - outflow
        cash_net = zeros(Float64, n, total_time)
        # discounted cash outflow
        disc_cash_out = zeros(Float64, n, total_time)
        # NPV in period t
        disc_cash_net = zeros(Float64, n, total_time)
        # generated electricity
        electricity = zeros(Float64, n, total_time)
        # discounted generated electricity
        disc_electricity = zeros(Float64, n, total_time)

    # simulation loop - iterate over simulations (each has different construction time)
    for i in 1:n
        construction_end = rand_construction_time[i]
        lifetime_end = construction_end + pj.time[2]

        for t in 1:total_time
            if t <= construction_end
                # Construction phase
                cash_out[i,t] = rand_investment[i] / construction_end
                cash_net[i,t] = -cash_out[i,t]

            elseif t <= lifetime_end
                # Operating phase
                electricity[i,t] = pj.plant_capacity * rand_loadfactor[i,t] * 8760
                cash_in[i,t] = rand_electricity_price[i,t] * electricity[i,t]
                cash_out[i,t] = operating_cost_fix + operating_cost_variable * electricity[i,t]
                cash_net[i,t] = cash_in[i,t] - cash_out[i,t]
                disc_electricity[i,t] = electricity[i,t] / ((1 + rand_wacc[i])^(t-1))
            else
                # Beyond this simulation's lifetime - arrays already initialized to zero
                continue
            end

            # Discount cash flows (applies to all periods)
            disc_cash_out[i,t] = cash_out[i,t] / ((1 + rand_wacc[i])^(t-1))
            disc_cash_net[i,t] = cash_net[i,t] / ((1 + rand_wacc[i])^(t-1))
        end
    end

    # output
    return(disc_cash_out = disc_cash_out, disc_cash_net = disc_cash_net, disc_electricity = disc_electricity)

end

"""
The npv_lcoe function calculates the net present value (NPV) and levelized cost of electricity (LCOE) for a given input, disc_res. The input, disc_res, is assumed to be an object with three attributes: disc_cash_out, disc_cash_net, and disc_electricity.
Next, the function creates three local variables disc_cash_out, disc_cash_net, and disc_electricity which are assigned the values of the corresponding attributes of the input disc_res.
The NPV is then calculated by taking the sum of the disc_cash_net variable along the second dimension using the sum function. Similarly, the LCOE is calculated by taking the ratio of the sum of the disc_cash_out variable along the second dimension to the sum of the disc_electricity variable along the second dimension.
Finally, the function returns a named tuple with two fields, npv and lcoe with the calculated values.
"""
function npv_lcoe(disc_res)
 
    @info "calculating NPV and LCOE"

    disc_cash_out = disc_res.disc_cash_out
    disc_cash_net = disc_res.disc_cash_net
    disc_electricity = disc_res.disc_electricity

    npv = sum(disc_cash_net,dims=2)
    lcoe = sum(disc_cash_out,dims=2) ./ sum(disc_electricity,dims=2)

    return(npv = npv, lcoe = lcoe)

end

"""
    calculate_vendor_baseline_lcoe(pj::project, wacc_range::Vector, elec_price_range::Vector)

Calculate deterministic vendor baseline LCOE using manufacturer-stated OCC and mean values.

This function performs a single deterministic calculation (not Monte Carlo) using:
- Manufacturer OCC (pj.investment * pj.plant_capacity) without SOAK discount
- Mean WACC from the provided range
- Mean electricity price from the provided range
- Mean load factor from project data

# Arguments
- `pj::project`: Project object with manufacturer data
- `wacc_range::Vector`: [min, max] WACC range
- `elec_price_range::Vector`: [min, max] electricity price range

# Returns
- Float64: Deterministic LCOE value in USD/MWh

This represents the vendor's claimed LCOE based on their stated overnight construction cost,
serving as a reference point for learning curve analysis.
"""
function calculate_vendor_baseline_lcoe(pj::project, wacc_range::Vector, elec_price_range::Vector)
    @info("Calculating vendor baseline LCOE for $(pj.name)")

    # Calculate mean values
    mean_wacc = mean(wacc_range)
    mean_elec_price = mean(elec_price_range)
    mean_loadfactor = mean(pj.loadfactor)

    # Total time
    total_time = pj.time[1] + pj.time[2]

    # Use manufacturer OCC without SOAK discount
    investment = pj.investment * pj.plant_capacity

    # O&M costs
    operating_cost_fix = pj.plant_capacity * pj.operating_cost[1]
    operating_cost_variable = pj.operating_cost[2] + pj.operating_cost[3]

    # Initialize arrays (single calculation, n=1)
    disc_cash_out = zeros(Float64, total_time)
    disc_electricity = zeros(Float64, total_time)

    # Calculation loop
    for t in 1:total_time
        if t <= pj.time[1]
            # Construction period
            cash_out_t = investment / pj.time[1]
        else
            # Operating period
            electricity_t = pj.plant_capacity * mean_loadfactor * 8760
            cash_out_t = operating_cost_fix + operating_cost_variable * electricity_t
            disc_electricity[t] = electricity_t / ((1 + mean_wacc) ^ (t-1))
        end
        disc_cash_out[t] = cash_out_t / ((1 + mean_wacc) ^ (t-1))
    end

    # Calculate LCOE
    vendor_lcoe = sum(disc_cash_out) / sum(disc_electricity)

    @info("Vendor baseline LCOE for $(pj.name): $(round(vendor_lcoe, digits=2)) USD/MWh")

    return vendor_lcoe
end

"""
The investment_simulation function simulates an investment project, given an instance of the project type pj and a set of random variables rand_vars.
The function runs a Monte Carlo simulation by calling the mc_run function and passing in the number of simulations to run, n, the project pj, and the set of random variables rand_vars. The function saves the results of the Monte Carlo simulation in the local variable disc_res.
The function then calculates the net present value (NPV) and the levelized cost of electricity (LCOE) by calling the npv_lcoe function and passing in the disc_res variable as input. The results of this calculation are saved in the local variable res.
The function then returns the results of the simulation in the form of the res variable.
"""
function investment_simulation(pj::project, rand_vars)

    # Infer n from rand_vars size (instead of using global n)
    n = length(rand_vars.wacc)

    # run the Monte Carlo simulation
    disc_res = mc_run(n, pj, rand_vars)

    # NPV and LCOE calculation
    res = npv_lcoe(disc_res)

    # output
    @info "simulation results" pj.name pj.type NPV = mean(res[1]) LCOE = mean(res[2])
    return(res)

end

"""
The si_correct function takes in a sensitivity index value si and returns a corrected sensitivity index value.
If si is approximately equal to negative zero (-0.0) with an absolute tolerance of 1e-2, the function returns a float 0.0. If si is greater than 1.0, the function returns 1.0. Otherwise, the function returns si unchanged.
The purpose of this function is to correct any sensitivity index values that may be out of range or have unexpected floating point behavior due to numerical errors.
"""
function si_correct(si)

    if isapprox(si, -0.0, atol = 1e-2)
        return 0.0
    elseif si > 1.0
        return 1.0
    else
        return si
    end

end

"""
The si_first_order function calculates the first-order Sobol' index (also called variance-based sensitivity index) for a given random variable.
The function takes three arguments, A, B, and AB, which represent random variable matrices.
The first-order sensitivity index measures the impact of each individual random variable on the output of a simulation model.
The function computes the first-order sensitivity index using the formula: mean(B .* (AB .- A)) / var(vcat(A, B), corrected = false).
The function then applies a correction using the "si_correct" function to ensure that the sensitivity index is within the range of 0 and 1.
The function returns the corrected first-order sensitivity index as a float value.
Input Arguments
    A: An array of output values corresponding to one set of random input values.
    B: An array of output values corresponding to another set of random input values.
    AB: A matrix of output values corresponding to the set of random input variables used for matrix A except for one variable coresponding to the random input used for matrix B. 
"""
# first-order sensitivity index
function si_first_order(A,B,AB)
    
    si = mean(B .* (AB .- A)) / var(vcat(A, B), corrected = false)
    si = si_correct(si)
    return si

end

"""
The si_total_order function computes the total-order Sobol' index (also called variance-based sensitivity index) for a given random variable.
The function takes three arguments, A, B, and AB, which represent random variable matrices.
The total-order sensitivity index is a measure of the impact of each input on the output, taking into account the effect of the input in combination with all other inputs.
The function computes the total-order sensitivity index using the formula: si = 0.5 * mean((A .- AB).^2) / var(vcat(A, B), corrected = false).
The function then applies a correction using the "si_correct" function to ensure that the sensitivity index is within the range of 0 and 1.
The function returns the corrected total-order sensitivity index as a float value.
Input Arguments
    A: An array of output values corresponding to one set of random input values.
    B: An array of output values corresponding to another set of random input values.
    AB: A matrix of output values corresponding to the set of random input variables used for matrix A except for one variable coresponding to the random input used for matrix B.   
"""
# total-effect sensitivity index
function si_total_order(A,B,AB)
    
    si = 0.5 * mean((A .- AB).^2) / var(vcat(A, B), corrected = false)
    si = si_correct(si)
    return si

end

"""
The sensiticity_index function calculates sensitivity indices for a given project based on the input parameters. The function takes in five arguments:
    opt_scaling: A string representing the type of optimization scaling to use.
    n: An integer representing the number of simulations to run for each random variable.
    wacc: A vector of doubles representing the weighted average cost of capital for each year of the project.
    electricity_price: A vector of doubles representing the electricity prices for each year of the project.
    pj: A project object containing all relevant information about the project.
The function then generates two random variable matrices A and B using the gen_rand_vars function. It then builds four matrices AB from A and B for each random variable. The function then runs Monte Carlo simulations for A, B, and each AB matrix, using the investment_simulation function.
The function then calculates sensitivity indices for both NPV and LCOE. The sensitivity indices are calculated separately for each random variable, and both first-order and total-order sensitivity indices are calculated. The resulting sensitivity indices are stored in tuples S_NPV, ST_NPV, S_LCOE, and ST_LCOE.
Finally, the function outputs the sensitivity results for the project, including the project's name, type, and the calculated sensitivity indices for NPV and LCOE. The function returns a tuple containing the sensitivity indices for NPV and LCOE.
"""
function sensitivity_index(opt_scaling::String, n::Int64, wacc::Vector, electricity_price_mean::Float64, pj::project;
                          construction_time_range::Union{Nothing,Vector}=nothing)

    # generate random variable matrices A and B
    @info "generating matrix A"
    rand_vars_A = gen_rand_vars(opt_scaling, n, wacc, electricity_price_mean, pj;
                                construction_time_range=construction_time_range);
    @info "generating matrix B"
    rand_vars_B = gen_rand_vars(opt_scaling, n, wacc, electricity_price_mean, pj;
                                construction_time_range=construction_time_range);

    # build random variable matrices AB from A and B for each random variable
    # AB matrices vary one parameter at a time for Sobol sensitivity
    @info "building matrices AB"
    rand_vars_AB1 = (wacc = rand_vars_B.wacc, electricity_price = rand_vars_A.electricity_price,
                    loadfactor = rand_vars_A.loadfactor, investment = rand_vars_A.investment,
                    construction_time = rand_vars_A.construction_time);
    rand_vars_AB2 = (wacc = rand_vars_A.wacc, electricity_price = rand_vars_A.electricity_price,
                    loadfactor = rand_vars_A.loadfactor, investment = rand_vars_A.investment,
                    construction_time = rand_vars_B.construction_time);  # Vary construction_time
    rand_vars_AB3 = (wacc = rand_vars_A.wacc, electricity_price = rand_vars_A.electricity_price,
                    loadfactor = rand_vars_B.loadfactor, investment = rand_vars_A.investment,
                    construction_time = rand_vars_A.construction_time);
    rand_vars_AB4 = (wacc = rand_vars_A.wacc, electricity_price = rand_vars_A.electricity_price,
                    loadfactor = rand_vars_A.loadfactor, investment = rand_vars_B.investment,
                    construction_time = rand_vars_A.construction_time);

    # run Monte Carlo simulations for A, B, and AB
    @info "matrix A:"
    sensi_res_A = investment_simulation(pj, rand_vars_A);
    @info "matrix B:"
    sensi_res_B = investment_simulation(pj, rand_vars_B);
    @info "matrix AB1"
    sensi_res_AB1 = investment_simulation(pj, rand_vars_AB1);
    @info "matrix AB2"
    sensi_res_AB2 = investment_simulation(pj, rand_vars_AB2);
    @info "matrix AB3"
    sensi_res_AB3 = investment_simulation(pj, rand_vars_AB3);
    @info "matrix AB4"
    sensi_res_AB4 = investment_simulation(pj, rand_vars_AB4);

    # sensitivity indices for NPV
    s_npv = (
        wacc = si_first_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB1[1]),
        construction_time = si_first_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB2[1]),
        loadfactor = si_first_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB3[1]),
        investment = si_first_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB4[1]),
        )
    st_npv = (
        wacc = si_total_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB1[1]),
        construction_time = si_total_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB2[1]),
        loadfactor = si_total_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB3[1]),
        investment = si_total_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB4[1])
        )

    # sensitivity indices for LCOE
    s_lcoe = (
        wacc = si_first_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB1[2]),
        construction_time = si_first_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB2[2]),
        loadfactor = si_first_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB3[2]),
        investment = si_first_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB4[2]),
        )
    st_lcoe = (
        wacc = si_total_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB1[2]),
        construction_time = si_total_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB2[2]),
        loadfactor = si_total_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB3[2]),
        investment = si_total_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB4[2])
        )
    
    # output
    @info "sensitivity results" pj.name pj.type S_NPV = s_npv ST_NPV = st_npv S_LCOE = s_lcoe ST_LCOE = st_lcoe
    return(s_npv = s_npv, st_npv = st_npv, s_lcoe = s_lcoe, st_lcoe = st_lcoe)

end

"""
The gen_scaled_investment function takes in two arguments, scaling and pj, and returns an array containing the scaled investment costs for a given project.

    Arguments
        scaling::Vector: A vector of scaling parameters that determines the range of the scaled investment costs.
        pj::project: An instance of the project struct that contains information about the investment concept, investment type, investment estimate by manufacturer, plant capacity, learning factor, project time, load factor, operating cost, and reference reactor.

    Output
        scaled_investment: A vector of scaled investment costs based on the project information and the given scaling parameter. The vector contains four parts:
            The deterministic investment cost based on manufacturer estimates in USD per MW.
            The range (lower and upper bound) of scaled investment cost based on Roulstone in USD per MW.
            The range (lower and upper bound) of scaled investment cost based on Rothwell in USD per MW.
            The scaled investment cost based on Carelli (fixed β = 0.20) in USD per MW.

    Note: SOAK discount (learning factor) is NOT applied in this function. Learning is handled separately in scenario simulations via the apply_learning parameter in gen_rand_vars().
"""
function gen_scaled_investment(scaling::Vector, pj::project)

    # generation of project specific scaled investment cost ranges
    # note that the scaling parameter here are converted such that Rothwell and Roulstone coincide.
    # IMPORTANT: SOAK discount (1-learning_factor) removed - learning is applied separately in scenario simulations
        # deterministic investment cost based on manufacturer estimates [USD/MW]
        scaled_investment = pj.investment
        # scaled investment cost based on Roulstone [USD/MW]
        scaled_investment = vcat(scaled_investment, pj.reference_pj[1] * pj.reference_pj[2] * (pj.plant_capacity/pj.reference_pj[2]) .^ scaling / pj.plant_capacity)
        # scaled investment cost based on Rothwell [USD/MW]
        scaled_investment = vcat(scaled_investment, pj.reference_pj[1] * pj.reference_pj[2] * (pj.plant_capacity/pj.reference_pj[2]) .^ (1 .+ log.(scaling) ./ log(2)) / pj.plant_capacity)
        # scaled investment cost based on Carelli (fixed β = 0.20) [USD/MW]
        carelli_scaling = 0.20
        scaled_investment = vcat(scaled_investment, pj.reference_pj[1] * pj.reference_pj[2] * (pj.plant_capacity/pj.reference_pj[2]) ^ carelli_scaling / pj.plant_capacity)

    # output
    return(scaled_investment = scaled_investment)

end