# SOAK Discount Removal - Learning Separation Summary

## Issue
The code was applying a 10% SOAK discount (`1 - learning_factor = 0.9`) in the `gen_scaled_investment()` plotting function, which contradicts the decision to separate learning completely from the base Monte Carlo simulation and apply it only in separate learning scenarios.

## Changes Made

### 1. Removed SOAK Discount from `gen_scaled_investment()`
**File**: `functions.jl` (lines 623, 625, 628)

**BEFORE** (Applied 10% SOAK discount):
```julia
# Line 623: Roulstone scaling
scaled_investment = vcat(scaled_investment,
    pj.reference_pj[1] * pj.reference_pj[2] * (1-pj.learning_factor) *
    (pj.plant_capacity/pj.reference_pj[2]) .^ scaling / pj.plant_capacity)

# Line 625: Rothwell scaling
scaled_investment = vcat(scaled_investment,
    pj.reference_pj[1] * pj.reference_pj[2] * (1-pj.learning_factor) *
    (pj.plant_capacity/pj.reference_pj[2]) .^ (1 .+ log.(scaling) ./ log(2)) / pj.plant_capacity)

# Line 628: Carelli scaling
scaled_investment = vcat(scaled_investment,
    pj.reference_pj[1] * pj.reference_pj[2] * (1-pj.learning_factor) *
    (pj.plant_capacity/pj.reference_pj[2]) ^ carelli_scaling / pj.plant_capacity)
```

**AFTER** (No SOAK discount applied):
```julia
# Line 623: Roulstone scaling
scaled_investment = vcat(scaled_investment,
    pj.reference_pj[1] * pj.reference_pj[2] *
    (pj.plant_capacity/pj.reference_pj[2]) .^ scaling / pj.plant_capacity)

# Line 625: Rothwell scaling
scaled_investment = vcat(scaled_investment,
    pj.reference_pj[1] * pj.reference_pj[2] *
    (pj.plant_capacity/pj.reference_pj[2]) .^ (1 .+ log.(scaling) ./ log(2)) / pj.plant_capacity)

# Line 628: Carelli scaling
scaled_investment = vcat(scaled_investment,
    pj.reference_pj[1] * pj.reference_pj[2] *
    (pj.plant_capacity/pj.reference_pj[2]) ^ carelli_scaling / pj.plant_capacity)
```

**Impact**: Investment comparison plot (`fig_investment_comparison_by_scale.pdf`) will now show ~11% higher investment costs for Roulstone, Rothwell, and Carelli methods (1.0 / 0.9 = 1.111).

### 2. Updated Function Documentation
**File**: `functions.jl`

- **Lines 602-617**: Updated `gen_scaled_investment()` docstring to clarify that SOAK discount is NOT applied and learning is handled separately
- **Lines 136-162**: Updated `gen_rand_vars()` docstring to clarify that learning is disabled by default and only applied in separate scenario simulations

## Verification - No Learning in Base Monte Carlo

### ✅ Base Monte Carlo Simulation (`run_simulation.jl`)
```julia
# Line 15
rand_vars = gen_rand_vars(opt_scaling, n, wacc, electricity_price_mean, pjs[p];
                          construction_time_range=construction_time_range)
```
- Does NOT pass `apply_soak_discount` → defaults to `false`
- Does NOT pass `apply_learning` → defaults to `false`
- **Result**: NO learning applied ✓

### ✅ Learning Scenarios (`smr-mcs-learning.jl`)
```julia
# Line 144
rand_vars = gen_rand_vars(opt_scaling, n, wacc, electricity_price_mean, pjs[p];
                          apply_learning=applyL,
                          N_unit=Nunit,
                          LR=LR,
                          kappa=kappa,
                          floor_m=floor_m,
                          construction_time_range=construction_time_range)
```
- Uses `apply_learning` parameter for controlled learning scenarios
- Learning rates: LR=5%, 10%, 15%
- Experience levels: N=1, 2, 4, 6, 12
- **Result**: Learning only in separate scenarios ✓

### ✅ Plotting Function (`gen_scaled_investment()`)
- **Before**: Hardcoded 10% SOAK discount in investment comparison plot
- **After**: No SOAK discount - matches base simulation approach
- **Result**: Consistent with separation of learning ✓

## Summary

| Component | SOAK Discount | Learning Curves | Status |
|-----------|---------------|-----------------|--------|
| Base Monte Carlo (`run_simulation.jl`) | ❌ Not applied | ❌ Not applied | ✅ Correct |
| Learning Scenarios (`smr-mcs-learning.jl`) | ❌ Not applied | ✅ Applied (controlled) | ✅ Correct |
| Investment Plot (`gen_scaled_investment()`) | ❌ Not applied (FIXED) | ❌ Not applied | ✅ Fixed |

## Expected Changes in Results

### Investment Comparison Plot
The investment estimates for Roulstone, Rothwell, and Carelli will increase by ~11%:
- **Old**: Investment × 0.9 (with SOAK discount)
- **New**: Investment × 1.0 (no SOAK discount)
- **Example**: A reactor with 5000 USD/kW → was 4500 USD/kW, now 5000 USD/kW

### Monte Carlo Simulation Results
- **No change** - base simulation was already correct (no learning applied)
- Learning scenarios remain unchanged (already using controlled learning parameters)

## Files Modified
1. `functions.jl`:
   - Removed `(1-pj.learning_factor)` from lines 623, 625, 628
   - Updated docstrings for `gen_scaled_investment()` and `gen_rand_vars()`

## Rationale
Learning is now **completely separated** from base analysis:
1. **Base simulations**: Show raw manufacturer/scaled estimates without learning
2. **Learning scenarios**: Apply controlled learning curves with explicit parameters
3. **Investment plots**: Consistent with base simulation (no learning)

This separation allows clearer comparison between FOAK costs (base) and learning scenarios.
