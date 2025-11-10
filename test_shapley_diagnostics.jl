using Pkg
Pkg.activate(pwd())
using Statistics

# Define paths required by data.jl
inputpath = "_input"
outputpath = "_output"

include("functions.jl")
include("data.jl")

# Test parameters
n_test = 500  # For total variance estimation
wacc = [0.04, 0.10]
ct_range = [3, 7]
scaling = [0.20, 0.75]

println("="^80)
println("SHAPLEY CONDITIONAL VARIANCE DIAGNOSTICS")
println("="^80)
println("\nThis test runs ONLY the validation diagnostics to debug V(S) estimation.")
println("\nParameters:")
println("  n_outer = 30 (outer loop iterations)")
println("  n_inner = 100 (inner loop iterations)")
println("  n_test = $n_test (for total variance estimation)")
println("="^80)

# Select test project
test_project = pjs[1]
electricity_price_mean = mean([52.2, 95.8])

println("\nTest Project: $(test_project.name)")
println("  Scale: $(test_project.scale)")
println("  Type: $(test_project.type)")
println("  Capacity: $(test_project.plant_capacity) MWe")

# Generate base samples for total variance estimation
println("\n" * "="^80)
println("Step 1: Estimating Total Variance")
println("="^80)

# Seed global RNG for reproducible A/B samples
using Random
Random.seed!(54321)
println("Global RNG seeded (54321) for reproducible A/B sample generation")

rand_vars_A = gen_rand_vars("rothwell", n_test, wacc, electricity_price_mean, test_project;
                            construction_time_range=ct_range)
rand_vars_B = gen_rand_vars("rothwell", n_test, wacc, electricity_price_mean, test_project;
                            construction_time_range=ct_range)

res_A = investment_simulation(test_project, rand_vars_A)
res_B = investment_simulation(test_project, rand_vars_B)

npv_all = vcat(res_A.npv, res_B.npv)
lcoe_all = vcat(res_A.lcoe, res_B.lcoe)
total_var_npv = var(npv_all, corrected=false)
total_var_lcoe = var(lcoe_all, corrected=false)

println("Total variance (from $(2*n_test) samples):")
println("  NPV:  $(round(total_var_npv, sigdigits=6))")
println("  LCOE: $(round(total_var_lcoe, sigdigits=6))")

# Diagnostic parameters
n_outer = 50
n_inner = 100

# SEED GLOBAL RNG for reproducibility
# This ensures outer_base is identical across runs
using Random
Random.seed!(12345)
println("\n⚠️  Global RNG seeded (12345) for reproducible outer_base generation")

# Generate shared outer design (same X_S^(k) for all coalitions)
println("Generating shared outer design (n_outer=$n_outer samples)...")
outer_base = []
for k in 1:n_outer
    sample = gen_rand_vars("rothwell", 1, wacc, electricity_price_mean, test_project;
                          construction_time_range=ct_range)
    push!(outer_base, sample)
end
println("Outer design generated: $(length(outer_base)) realizations")

println("\n" * "="^80)
println("CRITICAL DIAGNOSTICS - Run These to Find the Bug")
println("="^80)

# Test 1: V(∅) - Empty coalition
println("\n[Test 1] V(∅) - No parameters fixed")
println("Expected: < 1% of total variance")
println("Interpretation: If large, inner loop isn't varying parameters properly")
coalition_id_empty = hash(Symbol[])
V_empty_npv, V_empty_lcoe = estimate_conditional_variance_V_S(
    test_project, Symbol[], outer_base, n_inner,
    "rothwell", wacc, electricity_price_mean,
    ct_range, coalition_id_empty
)
println("  V(∅) NPV:  $(round(V_empty_npv, sigdigits=4))")
println("  Total var NPV: $(round(total_var_npv, sigdigits=4))")
println("  Ratio: $(round(V_empty_npv/total_var_npv, digits=4)) (MUST be < 0.01)")
if V_empty_npv/total_var_npv > 0.01
    println("  ❌ FAIL: V(∅) too large! Inner loop not working correctly.")
else
    println("  ✓ PASS: V(∅) is negligible")
end

println("\n  V(∅) LCOE: $(round(V_empty_lcoe, sigdigits=4))")
println("  Total var LCOE: $(round(total_var_lcoe, sigdigits=4))")
println("  Ratio: $(round(V_empty_lcoe/total_var_lcoe, digits=4)) (MUST be < 0.01)")
if V_empty_lcoe/total_var_lcoe > 0.01
    println("  ❌ FAIL: V(∅) too large! Inner loop not working correctly.")
else
    println("  ✓ PASS: V(∅) is negligible")
end

# Test 2: V(all) - All parameters fixed
println("\n" * "-"^80)
println("[Test 2] V({WACC, CT, LF, Inv}) - All parameters fixed")
println("Expected: 90-110% of total variance")
println("Interpretation: If too low, outer loop isn't fixing parameters properly")
coalition_id_full = hash([:wacc, :construction_time, :loadfactor, :investment])
V_all_npv, V_all_lcoe = estimate_conditional_variance_V_S(
    test_project, [:wacc, :construction_time, :loadfactor, :investment],
    outer_base, n_inner, "rothwell", wacc, electricity_price_mean,
    ct_range, coalition_id_full
)
println("  V(all) NPV:  $(round(V_all_npv, sigdigits=4))")
println("  Total var NPV: $(round(total_var_npv, sigdigits=4))")
println("  Ratio: $(round(V_all_npv/total_var_npv, digits=4)) (MUST be 0.9-1.1)")
if V_all_npv/total_var_npv < 0.9 || V_all_npv/total_var_npv > 1.1
    println("  ❌ FAIL: V(all) not close to total variance!")
    if V_all_npv/total_var_npv < 0.9
        println("         Outer loop isn't fixing parameters correctly")
    else
        println("         Variance is being amplified somehow")
    end
else
    println("  ✓ PASS: V(all) ≈ total variance")
end

println("\n  V(all) LCOE: $(round(V_all_lcoe, sigdigits=4))")
println("  Total var LCOE: $(round(total_var_lcoe, sigdigits=4))")
println("  Ratio: $(round(V_all_lcoe/total_var_lcoe, digits=4)) (MUST be 0.9-1.1)")
if V_all_lcoe/total_var_lcoe < 0.9 || V_all_lcoe/total_var_lcoe > 1.1
    println("  ❌ FAIL: V(all) not close to total variance!")
else
    println("  ✓ PASS: V(all) ≈ total variance")
end

# Test 3: Individual parameters
println("\n" * "-"^80)
println("[Test 3] V({param}) - Single parameter effects")
println("Expected LCOE contributions:")
println("  WACC: 30-45%, CT: 15-30%, LF: 5-15%, Inv: 10-25%")

V_individual_npv = Dict{Symbol, Float64}()
V_individual_lcoe = Dict{Symbol, Float64}()

for param in [:wacc, :construction_time, :loadfactor, :investment]
    println("\n  Computing V({$param})...")
    coalition_id_param = hash([param])
    V_p_npv, V_p_lcoe = estimate_conditional_variance_V_S(
        test_project, [param], outer_base, n_inner,
        "rothwell", wacc, electricity_price_mean,
        ct_range, coalition_id_param
    )
    V_individual_npv[param] = V_p_npv
    V_individual_lcoe[param] = V_p_lcoe

    pct_npv = round(100*V_p_npv/total_var_npv, digits=1)
    pct_lcoe = round(100*V_p_lcoe/total_var_lcoe, digits=1)
    println("  V({$param})")
    println("    NPV:  $(round(V_p_npv, sigdigits=4)) ($(pct_npv)% of total)")
    println("    LCOE: $(round(V_p_lcoe, sigdigits=4)) ($(pct_lcoe)% of total)")
end

# Test 4: Superadditivity check
println("\n" * "-"^80)
println("[Test 4] Superadditivity and Double-Counting Check")
sum_V_individual_npv = sum(values(V_individual_npv))
sum_V_individual_lcoe = sum(values(V_individual_lcoe))

println("\nFor INDEPENDENT parameters: Σ V({i}) ≤ V(all)")
println("For CORRELATED parameters: Σ V({i}) may exceed V(all) (double-counting)")
println("\nNPV:")
println("  Σ V({i}):     $(round(sum_V_individual_npv, sigdigits=4))")
println("  V(all):       $(round(V_all_npv, sigdigits=4))")
println("  Ratio:        $(round(sum_V_individual_npv/V_all_npv, digits=3))")
if sum_V_individual_npv > V_all_npv
    println("  → Double-counting detected (expected for correlated WACC×CT)")
else
    println("  → No double-counting (parameters are independent)")
end

println("\nLCOE:")
println("  Σ V({i}):     $(round(sum_V_individual_lcoe, sigdigits=4))")
println("  V(all):       $(round(V_all_lcoe, sigdigits=4))")
println("  Ratio:        $(round(sum_V_individual_lcoe/V_all_lcoe, digits=3))")
if sum_V_individual_lcoe > V_all_lcoe
    over_count = round(100*(sum_V_individual_lcoe/V_all_lcoe - 1), digits=1)
    println("  → Double-counting: $(over_count)% over V(all)")
    println("     This is WHY we need Shapley (distributes shared variance fairly)")
else
    println("  → No double-counting (parameters are independent)")
end

# Test 5: Monotonicity
println("\n" * "-"^80)
println("[Test 5] Monotonicity Check: V should increase with coalition size")
println("Computing V({WACC}) and V({WACC, CT})...")

coalition_id_wacc = hash([:wacc])
V_wacc_npv, V_wacc_lcoe = estimate_conditional_variance_V_S(
    test_project, [:wacc], outer_base, n_inner,
    "rothwell", wacc, electricity_price_mean,
    ct_range, coalition_id_wacc
)

coalition_id_wacc_ct = hash([:wacc, :construction_time])
V_wacc_ct_npv, V_wacc_ct_lcoe = estimate_conditional_variance_V_S(
    test_project, [:wacc, :construction_time], outer_base, n_inner,
    "rothwell", wacc, electricity_price_mean,
    ct_range, coalition_id_wacc_ct
)

println("\nNPV:")
println("  V({WACC}):     $(round(V_wacc_npv, sigdigits=4))")
println("  V({WACC,CT}):  $(round(V_wacc_ct_npv, sigdigits=4))")
println("  Difference:    $(round(V_wacc_ct_npv - V_wacc_npv, sigdigits=4))")
if V_wacc_ct_npv >= V_wacc_npv
    println("  ✓ PASS: Monotonic (adding CT increases conditional variance)")
else
    println("  ❌ FAIL: Non-monotonic! V decreased when adding parameter")
end

println("\nLCOE:")
println("  V({WACC}):     $(round(V_wacc_lcoe, sigdigits=4))")
println("  V({WACC,CT}):  $(round(V_wacc_ct_lcoe, sigdigits=4))")
println("  Difference:    $(round(V_wacc_ct_lcoe - V_wacc_lcoe, sigdigits=4))")
if V_wacc_ct_lcoe >= V_wacc_lcoe
    println("  ✓ PASS: Monotonic (adding CT increases conditional variance)")
else
    println("  ❌ FAIL: Non-monotonic! V decreased when adding parameter")
end

# Summary
println("\n" * "="^80)
println("DIAGNOSTIC SUMMARY")
println("="^80)

all_pass = true

println("\nTest 1 (V(∅) ≈ 0):")
if V_empty_npv/total_var_npv < 0.01 && V_empty_lcoe/total_var_lcoe < 0.01
    println("  ✓ PASS")
else
    println("  ❌ FAIL - Inner loop not varying parameters")
    all_pass = false
end

println("\nTest 2 (V(all) ≈ total_var):")
npv_ok = 0.9 <= V_all_npv/total_var_npv <= 1.1
lcoe_ok = 0.9 <= V_all_lcoe/total_var_lcoe <= 1.1
if npv_ok && lcoe_ok
    println("  ✓ PASS")
else
    println("  ❌ FAIL - Outer loop not fixing parameters correctly")
    all_pass = false
end

println("\nTest 3 (Individual V({param}) reasonable):")
println("  Check manually if percentages match expectations")

println("\nTest 4 (Superadditivity for correlated params):")
if sum_V_individual_lcoe > V_all_lcoe
    println("  ✓ Expected: Double-counting confirms WACC×CT correlation")
else
    println("  ⚠ Unexpected: No double-counting (copula may not be working)")
end

println("\nTest 5 (Monotonicity):")
if V_wacc_ct_npv >= V_wacc_npv && V_wacc_ct_lcoe >= V_wacc_lcoe
    println("  ✓ PASS")
else
    println("  ❌ FAIL - V decreased when adding parameters")
    all_pass = false
end

if all_pass
    println("\n" * "="^80)
    println("✓ ALL CRITICAL TESTS PASSED")
    println("Conditional variance estimation V(S) is working correctly!")
    println("Shapley implementation should produce valid results.")
    println("="^80)
else
    println("\n" * "="^80)
    println("❌ SOME TESTS FAILED")
    println("See failures above to diagnose the bug in V(S) estimation.")
    println("="^80)
end

println("\nDiagnostics complete. Expected runtime: ~5-10 minutes")
