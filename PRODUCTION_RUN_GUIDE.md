# Production Run Guide

This guide explains how to run the final production simulations with optimized parameters.

## Summary of Changes

### 1. Monte Carlo Simulation (Step 1)
- **Before**: 100k samples
- **After**: 1M samples (10× increase)
- **Reason**: Matches Uncertainties report methodology for publication-quality results
- **Runtime**: ~50 minutes for 41 reactors (was ~5 min)

### 2. Sensitivity Index (Step 2)
- **Status**: UNCHANGED (100k samples)
- **Reason**: Results are somewhat invalid due to parameter correlations (Sobol assumes independence)
- **Runtime**: ~10 minutes
- **Note**: Kept for comparison purposes, but Shapley is the main analysis

### 3. Shapley Analysis (Step 3)
- **Logging**: NOW FULLY QUIET
- **Runtime**: ~4-5 hours for 41 reactors
- **Fix**: Added quiet mode to helper function to prevent 2GB log files

## Why VS Code Was Crashing

**Root Cause**: Excessive logging creating 2GB+ log files

**The Problem**:
1. Main shapley function: 32+ coalition log messages × 41 reactors = 1,312+ messages ✓ FIXED
2. Helper function: 3 messages × 32 coalitions × 41 reactors = 3,936 messages ❌ WAS NOT FIXED
3. **Total**: ~5,000 log messages = 2GB file → VS Code crash

**The Solution**:
- Added `quiet` parameter to `estimate_conditional_variance_V_S()`
- Now only logs reactor-level progress (41 messages total)
- Log files stay <50MB instead of >2GB

## How to Run Production Simulations

### Option 1: Command Line (RECOMMENDED for Shapley)

Run each step from a terminal (Git Bash, PowerShell, or Julia REPL):

```bash
# Step 1: Monte Carlo (1M samples, ~50 min)
julia run_1_mcs.jl

# Step 2: Sensitivity Index (~10 min)
julia run_2_sensitivity.jl

# Step 3: Shapley (4-5 hours)
# Use dedicated script to avoid ANY VS Code issues:
./run_shapley_cli.sh
```

**Why use the CLI script for Shapley?**
- Runs in a separate process (not VS Code Julia kernel)
- Logs to timestamped file
- Can monitor progress: `tail -f shapley_run_*.log`
- Cannot crash VS Code
- Can run overnight/in background

### Option 2: Julia REPL (if you prefer)

```julia
# Step 1
include("run_1_mcs.jl")

# Step 2
include("run_2_sensitivity.jl")

# Step 3
include("run_3_shapley.jl")  # May still crash VS Code due to long runtime
```

### Option 3: VS Code (NOT recommended for Shapley)

Can work for Steps 1-2, but **Shapley may still crash** due to:
- Long runtime (4-5 hours)
- VS Code timeout limits
- Julia kernel memory accumulation

## Monitoring Progress

### For Command-Line Runs

```bash
# Watch real-time progress
tail -f shapley_run_*.log

# Check if process is still running
ps aux | grep julia

# Check recent log entries
tail -20 shapley_run_*.log
```

### Expected Output

```
[ Info: Running Shapley sensitivity for reactor 1/41: BWRX-300
[ Info: Shapley sensitivity results
│   pj.name = "BWRX-300"
│   Sh_NPV = (wacc = 0.034, construction_time = 0.0, loadfactor = 0.0, investment = 0.969)
│   Sh_LCOE = (wacc = 0.353, construction_time = 0.139, loadfactor = 0.0, investment = 0.504)
[ Info: Completed 1/41 reactors
...
[ Info: Completed 41/41 reactors
[ Info: Shapley sensitivity analysis complete for all reactors!
```

## Computational Costs

| Step | Samples | Runtime | Evaluations | Output |
|------|---------|---------|-------------|--------|
| MCS | 1M × 41 | ~50 min | 41M | mcs-*.csv |
| SI | 100k × 41 | ~10 min | ~4M | si-*.csv |
| Shapley | 100k base + 50×100 nested | ~4-5 hours | ~4.5M | shapley-*.csv |
| **Total** | | **~5-6 hours** | **~50M** | |

## Output Files

After completion, you'll have:

```
_output/
  mcs-npv_results-rothwell.csv          # 1M × 41 NPV samples
  mcs-lcoe_results-rothwell.csv         # 1M × 41 LCOE samples
  mcs-*_summary-rothwell.csv            # Summary statistics
  si-npv_results-rothwell.csv           # Sobol indices (NPV)
  si-lcoe_results-rothwell.csv          # Sobol indices (LCOE)
  shapley-npv_results-rothwell.csv      # Shapley effects (NPV)
  shapley-lcoe_results-rothwell.csv     # Shapley effects (LCOE)
```

## Validation Checks

After Shapley completes, verify:

1. **Efficiency property**: Sum of Shapley effects ≈ 1.0
2. **Non-negativity**: All effects ≥ 0 (or slightly negative due to numerical noise)
3. **Reasonable magnitudes**:
   - WACC: 0.3-0.5 (dominant)
   - Investment: 0.3-0.5 (dominant)
   - Construction Time: 0.1-0.2
   - Load Factor: 0.0-0.1

## Troubleshooting

### Shapley Still Crashes

Try these in order:

1. **Use CLI script** (already immune to VS Code):
   ```bash
   ./run_shapley_cli.sh
   ```

2. **Reduce n_outer** (if desperate):
   Edit `functions.jl:1227`:
   ```julia
   n_outer = 30  # Was 50, reduces runtime by 40%
   ```

3. **Process in batches**:
   Edit `run_3_shapley.jl:89` to process reactors in chunks:
   ```julia
   # Process reactors 1-10
   for p in 1:10  # Instead of eachindex(pjs)
   ```
   Run multiple times, changing the range each time.

### Memory Issues

If Julia runs out of memory:
```julia
# Add garbage collection between reactors
GC.gc()  # Add after line 116 in run_3_shapley.jl
```

### Log Files Too Large

If even with quiet mode logs are too large:
- Don't use `tee` in CLI script (writes to both screen and file)
- Just redirect: `julia run_3_shapley.jl > log.txt 2>&1`

## Final Notes

- **MCS at 1M samples**: Publication-ready precision for tail risk estimates
- **SI kept at 100k**: Invalid anyway due to correlations, kept for comparison
- **Shapley now quiet**: No more 2GB log files
- **CLI script recommended**: Most reliable way to run Shapley
- **Total runtime**: ~5-6 hours for complete analysis

Good luck with your production runs!
