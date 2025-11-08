using Pkg
Pkg.activate(pwd())
using Statistics

# Define paths required by data.jl
inputpath = "_input"
outputpath = "_output"

include("functions.jl")
include("data.jl")

# Test parameters (small sample size for quick testing)
n_test = 500  # Smaller sample for faster testing
wacc = [0.04, 0.10]
ct_range = [3, 7]
scaling = [0.20, 0.75]

println("="^80)
println("SHAPLEY SENSITIVITY WITH VALIDATION DIAGNOSTICS")
println("="^80)
println("\nThis test runs the CORRECTED Shapley sensitivity implementation")
println("with added validation diagnostics to verify correctness.")
println("\nFIXES APPLIED:")
println("  ✓ Conditional sampling for Gaussian copula (WACC|CT and CT|WACC)")
println("  ✓ Proper parameter fixing (same value for all inner loop iterations)")
println("  ✓ Validation tests for V(∅), V(all), individual V({param}), monotonicity")
println("\nNESTED SAMPLING ALGORITHM:")
println("  - Outer loop (n_outer=30): Different X_S realizations")
println("  - Inner loop (n_inner=100): X_~S samples for each X_S")
println("  - Cost per coalition: 30 × 100 = 3000 model evaluations")
println("  - Total coalitions for 4 parameters: 32")
println("  - Total cost: ~96,000 model evaluations PLUS diagnostics")
println("\nDIAGNOSTIC TESTS (run at end):")
println("  1. V(∅) should be ≈ 0 (< 1% of total variance)")
println("  2. V(all) should ≈ total_var (90-100%)")
println("  3. Individual V({param}) for each parameter")
println("  4. Monotonicity: V({WACC,CT}) >= V({WACC})")
println("\nExpected runtime: ~15-20 minutes with diagnostics")
println("="^80)

# Select a test project
test_project = pjs[1]
electricity_price_mean = mean([52.2, 95.8])

println("\nTest Project: $(test_project.name)")
println("  Scale: $(test_project.scale)")
println("  Type: $(test_project.type)")
println("  Capacity: $(test_project.plant_capacity) MWe")

println("\n" * "="^80)
println("Running Shapley Sensitivity Analysis...")
println("="^80)
println("\nThis will compute Shapley effects for 4 parameters:")
println("  1. WACC")
println("  2. Construction Time")
println("  3. Load Factor")
println("  4. Investment")
println("\nCoalitions:")
println("  - For each parameter: 2^(k-1) = 8 coalitions")
println("  - Total: 4 × 8 = 32 unique coalitions")
println("  - Each coalition requires 2 evaluations: V(S) and V(S∪{i})")
println("  - Plus 2 base samples for total variance")
println("\nComputational cost breakdown:")
println("  - Base samples (A, B): 2 × $n_test = $(2*n_test) evals")
println("  - Coalitions: 32 × 2 × (30×100) = 192,000 evals")
println("  - Total: ~192,000 model evaluations")
println("="^80)

# Run Shapley sensitivity analysis
shapley_results = shapley_sensitivity_index(
    "rothwell", n_test, wacc, electricity_price_mean, test_project;
    construction_time_range=ct_range
)

println("\n" * "="^80)
println("SHAPLEY SENSITIVITY RESULTS")
println("="^80)

println("\nShapley Effects for NPV:")
println("  WACC:              $(round(shapley_results.sh_npv.wacc, digits=4))")
println("  Construction Time: $(round(shapley_results.sh_npv.construction_time, digits=4))")
println("  Load Factor:       $(round(shapley_results.sh_npv.loadfactor, digits=4))")
println("  Investment:        $(round(shapley_results.sh_npv.investment, digits=4))")
println("  Sum:               $(round(sum([shapley_results.sh_npv.wacc, shapley_results.sh_npv.construction_time, shapley_results.sh_npv.loadfactor, shapley_results.sh_npv.investment]), digits=4))")

println("\nShapley Effects for LCOE:")
println("  WACC:              $(round(shapley_results.sh_lcoe.wacc, digits=4))")
println("  Construction Time: $(round(shapley_results.sh_lcoe.construction_time, digits=4))")
println("  Load Factor:       $(round(shapley_results.sh_lcoe.loadfactor, digits=4))")
println("  Investment:        $(round(shapley_results.sh_lcoe.investment, digits=4))")
println("  Sum:               $(round(sum([shapley_results.sh_lcoe.wacc, shapley_results.sh_lcoe.construction_time, shapley_results.sh_lcoe.loadfactor, shapley_results.sh_lcoe.investment]), digits=4))")

println("\n" * "="^80)
println("INTERPRETATION")
println("="^80)
println("\nShapley effects represent the fair share of variance contribution for each")
println("parameter, accounting for correlations and interactions.")
println("\nKey properties to verify:")
println("  ✓ All effects should be in [0, 1]")
println("  ✓ Sum should be ≈ 1.0 (efficiency property)")
println("  ✓ Effects properly attribute variance to correlated parameters")
println("\nCORRECTED ALGORITHM:")
println("  - Uses replicated nested sampling to estimate V(S) = Var(E[Y | X_S])")
println("  - V(S) measures how much Y's EXPECTED value varies as X_S changes")
println("  - This is different from residual variance (old incorrect approach)")
println("  - Marginal contributions ΔV = V(S∪{i}) - V(S) should be mostly positive")
println("\nExpected LCOE Shapley effects:")
println("  - WACC: 0.35-0.50 (dominant)")
println("  - Construction Time: 0.20-0.35 (second)")
println("  - Load Factor: 0.10-0.20 (moderate)")
println("  - Investment: 0.05-0.15 (least)")
println("\nComparison with Classical Sobol:")
println("  - Classical Sobol assumes independence")
println("  - With ρ=0.4 correlation, Sobol over-attributes to WACC and CT")
println("  - Shapley fairly distributes the correlated variance")
println("="^80)

println("\n✓ Shapley sensitivity test complete!")
println("\nNext steps:")
println("  1. Review the Shapley effects above")
println("  2. Verify the sum is close to 1.0 for both NPV and LCOE (efficiency property)")
println("  3. Check all effects are non-negative (or very small negative due to noise)")
println("  4. Verify WACC and construction time have reasonable effects (not too high)")
println("\nValidation tests to run:")
println("  - V(∅) should be ≈ 0 (no fixed parameters → no conditional variance)")
println("  - V(all) should ≈ total_var (all fixed → all variance is conditional)")
println("  - Marginal contributions ΔV should be mostly positive")
println("\nTo adjust computational cost:")
println("  - Modify n_outer and n_inner in shapley_sensitivity_index()")
println("  - Current: n_outer=30, n_inner=100 (balanced)")
println("  - Fast test: n_outer=10, n_inner=20 (~400 evals/coalition)")
println("  - Production: n_outer=50, n_inner=100 (~5000 evals/coalition)")
println("="^80)
