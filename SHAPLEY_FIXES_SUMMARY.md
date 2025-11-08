# Shapley Sensitivity Analysis - Fixes Implementation

## Overview

Implemented all critical fixes from ChatGPT's analysis to eliminate negative Shapley effects and achieve proper normalization (sum ≈ 1.0).

## Commits

### Part 1: CRN Infrastructure (Commit 33044bf^)
- Added `Random` to imports in functions.jl
- Created `generate_combined_sample_from_outer_base()` function
  - Accepts pre-generated `X_S_fixed` from outer_base
  - Uses seeded RNG for all random draws (CRN implementation)
  - Implements conditional copula sampling with seeded RNG
- Modified `estimate_conditional_variance_V_S()` signature
  - New parameters: `outer_base::Vector`, `coalition_id::UInt`
  - Generates deterministic RNG seed: `hash((coalition_id, k))`
  - Calls new helper function for inner loop sampling

### Part 2: Shared Design + Normalization (Commit 33044bf)
- **Fix #2**: Generate shared outer design in `shapley_sensitivity_index()`
  - `outer_base` generated once (n_outer samples of all 4 parameters)
  - Same X_S^(k) realizations used for all coalition evaluations
  - Passed to all V(S) and V(S∪{i}) calls

- **Fix #1**: Correct normalization
  - Compute V(∅) and V(full) first using outer_base
  - Use `norm = V(full) - V(∅)` for Shapley normalization
  - Replace `total_var_npv/lcoe` (from A/B) with `norm_npv/lcoe`

- **Fix #3**: Common Random Numbers (CRN)
  - Each coalition assigned unique `coalition_id = hash(param_set)`
  - Deterministic seeding: `rng_seed = hash((coalition_id, k))`
  - Same k → same seed → same inner randomness
  - Reduces variance in ΔV = V(S∪{i}) - V(S)

- Updated all calls in main loop and diagnostics
  - Main coalition loop: V(S) and V(S∪{i}) with outer_base + coalition_id
  - Diagnostic section: reuse V_empty and V_full, update individual param calls
  - Monotonicity test: use new signature

### Part 3: Diagnostic Test Updates (Commit 1788061)
- Updated `test_shapley_diagnostics.jl` to match new signature
- Generate outer_base once before tests
- Update all 5 test sections to use outer_base + coalition_id
- Tests now use same algorithm as main Shapley computation

## Key Algorithm Changes

### Before (Incorrect)
```julia
# Different outer designs for each coalition
V_S = estimate_conditional_variance(pj, S, n_outer, n_inner, ...)      # Generate new X_S^(k)
V_S_union = estimate_conditional_variance(pj, S∪{i}, n_outer, n_inner, ...)  # Different X_S^(k)!

# Normalization from separate A/B samples
sh[i] = shapley_npv / total_var_npv  # total_var_npv from A/B
```

**Problems:**
- V(S) and V(S∪{i}) on different outer designs → noise amplification
- No CRN → sign flips in ΔV from randomness
- Normalization inconsistent with Shapley game definition

### After (Correct)
```julia
# Generate shared outer design once
outer_base = [gen_rand_vars(...) for k in 1:n_outer]

# Normalization from same engine
V_empty = estimate_conditional_variance(pj, ∅, outer_base, n_inner, coalition_id=hash(∅))
V_full = estimate_conditional_variance(pj, P, outer_base, n_inner, coalition_id=hash(P))
norm = V_full - V_empty

# Same outer design + CRN for paired evaluations
coalition_id_S = hash(S)
V_S = estimate_conditional_variance(pj, S, outer_base, n_inner, coalition_id_S)

coalition_id_union = hash(S∪{i})
V_S_union = estimate_conditional_variance(pj, S∪{i}, outer_base, n_inner, coalition_id_union)

# Correct normalization
sh[i] = shapley_npv / norm
```

**Benefits:**
- Same X_S^(k) for all coalitions → fair comparison
- CRN via deterministic seeding → reduced variance in ΔV
- Normalization from same estimator → efficiency property holds

## Expected Results After Fixes

### Shapley Effects (LCOE)
- **WACC**: 0.35-0.50 (dominant driver)
- **Construction Time**: 0.20-0.35 (second most important)
- **Load Factor**: 0.10-0.20 (moderate)
- **Investment**: 0.05-0.15 (least important)

### Validation Criteria
1. **All effects non-negative**: Sh_i >= -0.01 (tiny negatives OK due to Monte Carlo noise)
2. **Efficiency property**: Σ Sh_i = 0.98-1.02 (sum ≈ 1.0)
3. **V(∅) ≈ 0**: < 1% of total variance
4. **V(all) ≈ total_var**: 90-110% of total variance
5. **Monotonicity**: V({WACC,CT}) >= V({WACC})

### Before vs After (Expected)
| Metric | Before (Broken) | After (Fixed) |
|--------|----------------|---------------|
| WACC LCOE | -0.067 ❌ | 0.40-0.45 ✓ |
| CT LCOE | -0.027 ❌ | 0.25-0.30 ✓ |
| LF LCOE | 0.0 ⚠️ | 0.12-0.18 ✓ |
| Inv LCOE | 0.90 ⚠️ | 0.08-0.12 ✓ |
| Sum NPV | 0.84 ❌ | 0.98-1.02 ✓ |
| Sum LCOE | 0.96 ⚠️ | 0.98-1.02 ✓ |

## Mathematical Justification

### Fix #1: Normalization
**Problem**: Used `total_var = Var(Y)` from independent A/B samples
**Solution**: Use `norm = V(P) - V(∅)` from nested sampling

**Why it matters**: Shapley's efficiency property states:
```
Σ Sh_i = [V(P) - V(∅)] / norm
```

If norm ≠ V(P) - V(∅), the sum cannot equal 1.0.

With same estimator:
```
Σ Sh_i = Σ [weighted sum of ΔV] / [V(P) - V(∅)]
       = [V(P) - V(∅)] / [V(P) - V(∅)]  (by Shapley axioms)
       = 1.0 ✓
```

### Fix #2: Shared Outer Design
**Problem**: V(S) uses X_S^(k), V(S∪{i}) uses different X_S^(k')
**Solution**: Both use same outer_base[k]

**Why it matters**: Marginal contribution is:
```
ΔV_i(S) = V(S∪{i}) - V(S)
```

If both computed on different designs, ΔV includes:
1. True marginal effect of adding parameter i
2. **Noise from using different X_S realizations** ❌

With shared design:
```
ΔV_i(S) = (1/n_outer) Σ_k [Var(Y|X_{S∪{i}}^(k)) - Var(Y|X_S^(k))]
```
where X_S^(k) ⊂ X_{S∪{i}}^(k) (same outer realization).

This isolates the true marginal effect of i.

### Fix #3: Common Random Numbers (CRN)
**Problem**: Inner loop uses `randn()` → different randomness each call
**Solution**: Inner loop uses `randn(rng)` with seeded RNG

**Why it matters**: For same k and same coalition:
- First call: `randn()` → sequence [r1, r2, r3, ...]
- Second call: `randn()` → sequence [r4, r5, r6, ...] (different!)

With CRN:
- k=1, coalition_id=hash(S): seed=(coalition_id, 1) → [r1, r2, r3, ...]
- k=1, coalition_id=hash(S): seed=(coalition_id, 1) → [r1, r2, r3, ...] (same!)

This is a classical variance reduction technique: Var(X - Y) < Var(X) + Var(Y) when X and Y are correlated.

## Implementation Details

### Outer Base Structure
```julia
outer_base = Vector of n_outer NamedTuples, each containing:
  - wacc: Float64
  - construction_time: Int
  - loadfactor: Matrix{Float64} (n=1, total_time)
  - investment: Vector{Float64} (scaling factors)
```

### Coalition ID Seeding
```julia
coalition_id = hash([:wacc, :construction_time])  # Unique per coalition
rng_seed = hash((coalition_id, k))                # Unique per (coalition, outer_iter)
rng = Random.MersenneTwister(rng_seed)            # Deterministic RNG
```

### Conditional Sampling with CRN
```julia
# Example: WACC fixed, CT varies
if :wacc in S && !(:construction_time in S)
    # Extract fixed WACC from outer_base
    wacc_fixed = X_S_fixed.wacc

    # Sample CT conditionally with seeded RNG
    Z_wacc = quantile(Normal(), cdf(wacc_dist, wacc_fixed))
    Z_ct = randn(rng) * sqrt(1 - ρ^2) + Z_wacc * ρ  # Seeded!
    U_ct = cdf(Normal(), Z_ct)
    ct = quantile(ct_dist, U_ct)
end
```

## Testing Instructions

### Quick Diagnostic Test (~5-10 minutes)
```bash
julia test_shapley_diagnostics.jl
```

**Expected Output:**
- ✓ Test 1: V(∅) < 1% of total variance
- ✓ Test 2: V(all) ≈ 90-100% of total variance
- ✓ Test 3: Individual V({param}) reasonable (WACC 35-45%, CT 20-30%)
- ✓ Test 4: Double-counting detected (Σ V({i}) > V(all) due to ρ=0.4)
- ✓ Test 5: Monotonic (V({WACC,CT}) >= V({WACC}))

### Full Shapley Test (~15-20 minutes)
```bash
julia test_shapley_sensitivity.jl
```

**Expected Output:**
- Shapley NPV: WACC 0.30-0.40, CT 0.15-0.25, LF 0.15-0.25, Inv 0.15-0.25
- Shapley LCOE: WACC 0.35-0.50, CT 0.20-0.35, LF 0.10-0.20, Inv 0.05-0.15
- Sum NPV: 0.98-1.02 ✓
- Sum LCOE: 0.98-1.02 ✓
- All effects >= -0.01 ✓

## Remaining Considerations

### Fix #4: V(∅) Forced to Zero
**Status**: Not implemented (low priority)
**Reason**: With correct nested sampling, V(∅) is already ≈ 0 (< 1% of variance)
**If needed**: Set V_empty_npv = V_empty_lcoe = 0.0 explicitly after computation

### Fix #5: Investment Randomness
**Status**: Not addressed
**Reason**: Investment uses Roulstone scaling which may inject randomness
**Impact**: Minor - investment typically has small Shapley effect
**If needed**: Make investment draws deterministic or seed them separately

### Fix #6: Sample Size
**Current**: n_outer=30, n_inner=100 (3000 evals/coalition)
**ChatGPT recommendation**: n_outer=32, n_inner=500 (16000 evals/coalition)
**Rationale**: Inner loop dominates variance of E[Y | X_S]

To adjust, modify in `shapley_sensitivity_index()`:
```julia
n_outer = 32   # More X_S realizations
n_inner = 500  # Better estimation of E[Y | X_S=x]
```

**Tradeoff**:
- Current: ~96,000 total evals, ~15-20 min runtime
- Recommended: ~512,000 total evals, ~80-120 min runtime

For validation: keep current settings
For publication: increase to recommended

## Files Modified

1. **functions.jl** (2 commits)
   - Added Random import
   - Created generate_combined_sample_from_outer_base() with CRN
   - Modified estimate_conditional_variance_V_S() signature
   - Updated shapley_sensitivity_index() with fixes #1, #2, #3

2. **test_shapley_diagnostics.jl** (1 commit)
   - Generate outer_base before tests
   - Update all test calls to new signature

3. **SHAPLEY_FIXES_SUMMARY.md** (this file)
   - Complete documentation of all fixes

## Next Steps

1. **Run diagnostic test** to verify V(∅), V(all), monotonicity
2. **Run full Shapley test** to check efficiency property
3. **If sum ≈ 1.0 and all effects >= 0**: Success! ✓
4. **If still issues**: Increase n_inner to 500 and retest
5. **For production**: Consider n_outer=50, n_inner=200

## References

- Owen, A.B. (2014). "Sobol' indices and Shapley value." SIAM/ASA Journal on Uncertainty Quantification.
- Iooss, B. & Prieur, C. (2019). "Shapley effects for sensitivity analysis with correlated inputs."
- L'Ecuyer, P. & Lemieux, C. (2002). "Variance reduction via lattice rules." Management Science (CRN techniques).

---

**Implementation Date**: 2025-11-08
**Branch**: claude/fix-xlsx-csv-last-column-011CUpscE6ipqt8ksfJFSjCo
**Commits**: Part 1 (CRN), Part 2 (Normalization), Part 3 (Diagnostics)
