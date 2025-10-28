# Learning Curves Implementation

This document describes the learning curve functionality for modeling cost reductions through experience (learning-by-doing) in reactor manufacturing.

## Overview

The learning curve model implements the standard power law relationship where costs decrease as cumulative production increases:

```
Cost(N) = κ × SOAK_Cost × (1 - LR)^(log₂(N))
```

Where:
- **N**: Number of units built (experience level)
- **LR**: Learning rate (fraction cost reduction per doubling, e.g., 0.10 = 10%)
- **κ (kappa)**: FOAK premium multiplier (e.g., 1.20 = 20% premium above SOAK)
- **SOAK_Cost**: Standard-Of-A-Kind baseline cost (after initial learning)

## Key Components

### 1. Helper Function: `learning_multiplier()`

Located in: `functions.jl`

```julia
learning_multiplier(N::Int, LR::Float64; kappa::Float64=1.0, floor::Union{Nothing,Float64}=nothing)
```

**Parameters:**
- `N`: Number of units built
- `LR`: Learning rate (e.g., 0.10 for 10% per doubling)
- `kappa`: FOAK premium (default 1.0 = no premium)
- `floor`: Optional minimum multiplier (e.g., 1.0 to stop at SOAK)

**Returns:** Multiplier to apply to investment cost

**Examples:**
```julia
# FOAK with 20% premium, LR=10%
learning_multiplier(1, 0.10, kappa=1.20)  # → 1.20 (20% above SOAK)

# Second unit
learning_multiplier(2, 0.10, kappa=1.20)  # → 1.08 (~8% above SOAK)

# Fourth unit with floor at SOAK
learning_multiplier(4, 0.10, kappa=1.20, floor=1.0)  # → 1.00 (hits SOAK floor)

# Eighth unit with floor
learning_multiplier(8, 0.10, kappa=1.20, floor=1.0)  # → 1.00 (stays at SOAK)
```

### 2. Modified Function: `gen_rand_vars()`

The investment random variable generator now accepts optional learning parameters:

```julia
gen_rand_vars(opt_scaling, n, wacc, electricity_price, pj;
              apply_learning=false,      # Enable learning curve
              N_unit=1,                  # Number of units built
              LR=0.0,                    # Learning rate
              kappa=1.0,                 # FOAK premium multiplier
              floor_m=nothing)           # Optional floor
```

**Default behavior:** When `apply_learning=false`, the function behaves exactly as before (backward compatible).

**Learning application:**
- Applied AFTER all existing scaling logic (Manufacturer/Roulstone/Rothwell)
- Excluded for Large reactors by default (only applies to SMR/Micro)
- Multiplies `rand_investment` by the learning multiplier

### 3. Learning Scenario Runner: `smr-mcs-learning.jl`

This script runs multiple learning scenarios in one execution:

**Default scenarios (simplified):**
1. **baseline**: No learning applied (reference case)
2. **FOAK**: First-Of-A-Kind with 20% premium (N=1, LR=10%, κ=1.20)
3. **SOAK**: Standard-Of-A-Kind after learning (N=4, reaches baseline)

**Output file count:** 6 files total (3 scenarios × 2 files each)

**Usage:**
```bash
julia smr-mcs-learning.jl
```

**Configuration options in script:**
```julia
# Set to true to also save NPV files (creates 12 files instead of 6)
save_npv_files = false  # Default: only LCOE files

# Uncomment extended scenarios for full learning curve (creates 24 files)
# learning_cases = [baseline, N1, N2, N4, N6, N8]
```

**Output files (default):**
- `mcs-lcoe_results-{scaling}-{tag}.csv`: Full Monte Carlo LCOE results
- `mcs-lcoe_summary-{scaling}-{tag}.csv`: LCOE summary statistics

**Additional files (if save_npv_files=true):**
- `mcs-npv_results-{scaling}-{tag}.csv`: NPV results
- `mcs-npv_summary-{scaling}-{tag}.csv`: NPV summary

### 4. Plotting Functions: `functions_plots.jl`

Two new plotting functions:

#### a) `learning_curve_plot()`

Plot mean LCOE vs. N for learning scenarios:

```julia
# Default scenarios (3 cases)
learning_scenarios = [(1, "baseline"), (1, "FOAK"), (4, "SOAK")]

# Or extended scenarios (if you uncommented them in smr-mcs-learning.jl)
# learning_scenarios = [(1, "baseline"), (1, "LR10_N1_k120"), (2, "LR10_N2_k120"),
#                       (4, "LR10_N4_k120"), (6, "LR10_N6_k120"), (8, "LR10_N8_k120")]

fig = learning_curve_plot("_output", "roulstone", learning_scenarios)
save("_output/fig-learning_curve.pdf", fig)
```

**Options:**
- `reactor_name`: Plot specific reactor (default: average across all)
- `scale_filter`: Filter by scale ("Micro", "SMR", "Large")

#### b) `learning_curve_comparison_plot()`

Multi-panel comparison across reactor scales:

```julia
fig = learning_curve_comparison_plot("_output", "roulstone", learning_scenarios, pjs)
save("_output/fig-learning_curve_comparison.pdf", fig)
```

Shows three subplots (Micro, SMR, Large) with separate learning curves.

## Example Workflow

### Step 1: Run Learning Scenarios

```bash
julia smr-mcs-learning.jl
```

This generates CSV files for each scenario in `_output/`.

### Step 2: Create Plots

```julia
# Load the plotting functions
include("functions_plots.jl")
include("data.jl")  # For pjs vector

# Define scenarios to plot (matching what you ran in smr-mcs-learning.jl)
learning_scenarios = [
    (1, "baseline"),
    (1, "FOAK"),
    (4, "SOAK")
]

# Overall learning curve
fig1 = learning_curve_plot("_output", "roulstone", learning_scenarios)
save("_output/fig-learning_curve.pdf", fig1)

# By scale comparison
fig2 = learning_curve_comparison_plot("_output", "roulstone", learning_scenarios, pjs)
save("_output/fig-learning_curve_by_scale.pdf", fig2)
```

## Unit Test Examples

### Test 1: Basic Multiplier Calculation

```julia
using Test

@testset "Learning Multiplier Tests" begin
    # FOAK premium (N=1)
    @test learning_multiplier(1, 0.10, kappa=1.20) ≈ 1.20

    # Second unit (one doubling)
    @test learning_multiplier(2, 0.10, kappa=1.20) ≈ 1.08

    # Fourth unit (two doublings)
    @test learning_multiplier(4, 0.10, kappa=1.20) ≈ 0.972

    # Floor at SOAK (1.0)
    @test learning_multiplier(4, 0.10, kappa=1.20, floor=1.0) == 1.0
    @test learning_multiplier(8, 0.10, kappa=1.20, floor=1.0) == 1.0
end
```

### Test 2: Verify Capital-Only Effect

Run a test case and verify:
- Construction period costs change ✓
- Operating period costs unchanged ✓
- Only `rand_investment` is scaled

```julia
# Run baseline
rand_vars_baseline = gen_rand_vars("manufacturer", 100, wacc, elec_price, pj)

# Run with learning
rand_vars_learning = gen_rand_vars("manufacturer", 100, wacc, elec_price, pj;
                                   apply_learning=true, N_unit=2, LR=0.10, kappa=1.20)

# Check that investment changed
@test mean(rand_vars_learning.investment) < mean(rand_vars_baseline.investment)

# Check that other variables unchanged
@test rand_vars_learning.wacc == rand_vars_baseline.wacc
@test rand_vars_learning.electricity_price == rand_vars_baseline.electricity_price
@test rand_vars_learning.loadfactor == rand_vars_baseline.loadfactor
```

## Rationale for Design Choices

### 1. Why exclude Large reactors by default?

Learning-by-duplication primarily applies to:
- **Factory manufacturing** (SMR modules)
- **Serial production** (Micro reactors)
- **Standardized designs**

Large reactors are typically:
- **Site-built** (not factory-produced)
- **Custom designs** (limited standardization)
- **Low volume** (few units per site)

Therefore, learning curves are more applicable to SMR/Micro scales.

### 2. Why use a floor parameter?

The floor prevents unrealistic cost reductions:
- **floor=1.0**: Learning only burns off FOAK premium, stops at SOAK
- **floor=nothing**: Learning continues to NOAK (Nth-Of-A-Kind) levels
- **Realistic range**: Most studies suggest floors between 0.7-1.0

### 3. Why apply learning AFTER scaling?

Design rationale:
1. Scaling methods (Roulstone/Rothwell) capture **size effects**
2. Learning curves capture **experience effects**
3. These are independent phenomena
4. Applying sequentially: `Cost = BaseScale × SizeEffect × LearningEffect`

## Parameter Recommendations

Based on nuclear industry literature:

| Parameter | Typical Range | Conservative | Optimistic |
|-----------|---------------|--------------|------------|
| **LR** (Learning Rate) | 5-15% | 5% | 15% |
| **κ (FOAK premium)** | 1.15-1.30 | 1.30 | 1.15 |
| **Floor** | 0.7-1.0 | 1.0 (SOAK) | 0.7 (deep NOAK) |
| **N (units)** | 1-10 | 1-4 | 6-10 |

**Example scenarios:**

1. **Conservative**: LR=5%, κ=1.30, floor=1.0, N=1-4
2. **Base case**: LR=10%, κ=1.20, floor=1.0, N=1-6
3. **Optimistic**: LR=15%, κ=1.15, floor=0.8, N=1-10

## References

- Rothwell, G. (2015). *Economics of Nuclear Power*
- Grubler, A. (2010). "The costs of the French nuclear scale-up"
- Rubin et al. (2015). "A review of learning rates for electricity supply technologies"
- OECD/NEA (2020). "Unlocking Reductions in the Construction Costs of Nuclear"

## Troubleshooting

**Issue:** Learning not applied to Large reactors
- **Solution:** This is intentional. Remove the `pj.scale != "Large"` condition in `gen_rand_vars()` if needed.

**Issue:** Negative multipliers or unexpected values
- **Check:** LR should be between 0.0 and 0.3 (0-30%)
- **Check:** κ should be ≥ 1.0 (cannot have FOAK discount)
- **Check:** floor should be ≤ 1.0 (cannot floor above SOAK)

**Issue:** Plots show no learning effect
- **Check:** Verify `apply_learning=true` in scenario definitions
- **Check:** Confirm output files exist for all scenarios
- **Check:** Ensure N_unit values are sequential (1, 2, 4, 6, 8...)
