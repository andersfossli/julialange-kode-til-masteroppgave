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
println("SHAPLEY SENSITIVITY TEST")
println("="^80)
println("\nThis test verifies the Shapley sensitivity implementation works correctly")
println("with correlated inputs (WACC × construction_time, ρ=0.4)")
println("\nNote: This uses a small sample size (n=$n_test) for quick testing.")
println("Production runs should use n≥10000 for accurate results.")
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
println("\nFor 4 parameters, this requires 2^4 = 16 subset evaluations per parameter.")
println("Total: 4 × 16 = 64 model evaluations (plus base samples)")
println("\nExpected runtime: ~1-2 minutes with n=$n_test")
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
println("\nComparison with Classical Sobol:")
println("  - Classical Sobol indices assume independence")
println("  - With ρ=0.4 correlation between WACC and construction_time,")
println("    classical Sobol would over-attribute variance to these parameters")
println("  - Shapley effects fairly distribute the correlated variance")
println("="^80)

println("\n✓ Shapley sensitivity test complete!")
println("\nNext steps:")
println("  1. Review the Shapley effects above")
println("  2. Verify the sum is close to 1.0 for both NPV and LCOE")
println("  3. If tests pass, use shapley_sensitivity_index() for production runs")
println("  4. For accurate results, use n≥10000 in production")
println("="^80)
