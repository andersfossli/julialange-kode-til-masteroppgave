# Split Workflow for Production Runs

This workflow splits the full simulation into three independent steps to:
- Prevent VS Code crashes from massive log files
- Allow reviewing results at each stage before continuing
- Save time if later steps fail

## Why Split?

The original `smr-mcs.jl` creates **2GB+ log files** due to verbose Shapley logging (1,312+ coalition messages). This crashes VS Code and makes debugging impossible.

## Three-Step Workflow

### Step 1: Monte Carlo Simulation (~1-2 hours)

```bash
julia run_1_mcs.jl
```

**Runtime:** ~1-2 hours with 100k samples
**Output:**
- `mcs-npv_results-rothwell.csv` (41 reactors × 100k samples)
- `mcs-lcoe_results-rothwell.csv`
- `mcs-*_summary-rothwell.csv`
- `mcs-wacc_values-rothwell.csv`
- `mcs-investment_values-rothwell.csv`

**What it does:**
- Generates 100k random samples for each reactor
- Calculates NPV and LCOE distributions
- Saves raw results and summary statistics

### Step 2: Sobol Sensitivity Index (~30-60 minutes)

```bash
julia run_2_sensitivity.jl
```

**Runtime:** ~30-60 minutes
**Output:**
- `si-npv_results-rothwell.csv`
- `si-lcoe_results-rothwell.csv`

**What it does:**
- Computes first-order (S) and total-order (ST) Sobol indices
- Shows variance contribution of each parameter
- Fast variance-based sensitivity analysis

### Step 3: Shapley Sensitivity Analysis (~13-14 hours)

```bash
julia run_3_shapley.jl
```

**Runtime:** ~13-14 hours (41 reactors × 20 min each)
**Output:**
- `shapley-npv_results-rothwell.csv`
- `shapley-lcoe_results-rothwell.csv`

**What it does:**
- Computes Shapley effects for correlated parameters
- Uses **QUIET MODE** to suppress coalition-level logging
- Prevents 2GB log file issue

**Features:**
- Quiet mode: Only logs reactor-level progress (not 1,312 coalitions)
- Can be interrupted and resumed (just comment out completed reactors in the loop)
- Produces manageable log files (<100MB)

## Parameters

All scripts use the same production parameters:
- **n = 100,000** samples (10× more robust than 10k)
- **n_outer = 50** (Shapley outer loop)
- **n_inner = 100** (Shapley inner loop)
- **scaling = "rothwell"**

## Comparison: Original vs Split

| Aspect | Original (`smr-mcs.jl`) | Split Workflow |
|--------|------------------------|----------------|
| Runtime | ~24 hours total | Same (~24 hours) |
| Log file size | **2GB+** (crashes VS Code) | <100MB per step |
| Debugging | Impossible (too verbose) | Easy (step-by-step) |
| Resumability | Must restart from scratch | Resume at any step |
| Review | Wait 24 hours for all results | Review after each step |

## Quick Start

Run all three steps sequentially:

```bash
# Step 1: Monte Carlo (~1-2 hours)
julia run_1_mcs.jl

# Review MCS results before continuing
# Check: mcs-lcoe_results-rothwell.csv

# Step 2: Sensitivity Index (~30-60 minutes)
julia run_2_sensitivity.jl

# Review SI results before continuing
# Check: si-lcoe_results-rothwell.csv

# Step 3: Shapley (~13-14 hours)
julia run_3_shapley.jl

# Review final Shapley results
# Check: shapley-lcoe_results-rothwell.csv
```

## Troubleshooting

**If Step 3 (Shapley) takes too long:**
- Reduce `n` from 100k to 10k (10× faster, but less robust)
- Reduce `n_outer` from 50 to 30 (40% faster, but higher variance)
- Run overnight and check progress in the morning

**If you want verbose logging:**
Edit `run_3_shapley.jl` line 77:
```julia
quiet=true  # Change to quiet=false for full coalition logging
```
**Warning:** This will create a 2GB+ log file again!

**If Julia runs out of memory:**
- Process reactors in batches (split the `for p in eachindex(pjs)` loop)
- Save intermediate results after each batch

## Verification

After all steps complete, verify:
1. **MCS**: Check summary statistics are reasonable
2. **SI**: Check S + ST ≈ 1.0 for dominant parameters
3. **Shapley**: Check sum of effects ≈ 1.0 (efficiency property)
