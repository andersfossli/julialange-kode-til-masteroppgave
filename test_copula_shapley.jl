using Pkg
Pkg.activate(pwd())
using Statistics, Plots

# Define paths required by data.jl
inputpath = "_input"
outputpath = "_output"

include("functions.jl")
include("data.jl")

# Test parameters
n_test = 1000
wacc = [0.04, 0.10]
ct_range = [3, 7]

# Scaling parameter (required by rothwell/roulstone scaling methods)
scaling = [0.20, 0.75]

println("="^80)
println("TEST 1: Gaussian Copula with Triangular Marginals")
println("="^80)

# Generate samples
wacc_samples, ct_samples = generate_correlated_samples(n_test, wacc, ct_range; ρ=0.4)

# Verify correlation
actual_corr = cor(wacc_samples, ct_samples)
println("Target ρ: 0.4")
println("Actual ρ: $(round(actual_corr, digits=4))")
println("Deviation: $(round(abs(actual_corr - 0.4), digits=4))")

# Verify marginals
println("\nWACC distribution:")
println("  Range: [$(wacc[1]), $(wacc[2])]")
println("  Sample: [$(round(minimum(wacc_samples), digits=4)), $(round(maximum(wacc_samples), digits=4))]")
println("  Mean: $(round(mean(wacc_samples), digits=4)) (expected ≈ $(round(mean(wacc), digits=4)))")

println("\nConstruction Time distribution:")
println("  Range: [$(ct_range[1]), $(ct_range[2])]")
println("  Sample: [$(round(minimum(ct_samples), digits=4)), $(round(maximum(ct_samples), digits=4))]")
println("  Mean: $(round(mean(ct_samples), digits=4)) (expected ≈ $(round(mean(ct_range), digits=4)))")

# Plot
scatter(wacc_samples, ct_samples, alpha=0.3,
        xlabel="WACC", ylabel="Construction Time",
        title="Copula Test (ρ=$(round(actual_corr, digits=3)))",
        legend=false)
savefig("_output/test_copula.png")

println("\n✓ Copula test complete. Plot saved to _output/test_copula.png")

println("\n"*"="^80)
println("TEST 2: Triangular Distribution Validation")
println("="^80)

# Test triangular distribution properties
using Distributions
lf_min = 0.65
lf_max = 0.95
lf_mode = lf_min + 0.6 * (lf_max - lf_min)
lf_dist = TriangularDist(lf_min, lf_max, lf_mode)
lf_samples = rand(lf_dist, n_test)

println("\nLoad Factor Triangular Distribution:")
println("  Range: [$lf_min, $lf_max]")
println("  Mode: $lf_mode")
println("  Sample: [$(round(minimum(lf_samples), digits=4)), $(round(maximum(lf_samples), digits=4))]")
println("  Mean: $(round(mean(lf_samples), digits=4)) (expected ≈ $(round(mean(lf_dist), digits=4)))")
println("  Std Dev: $(round(std(lf_samples), digits=4)) (expected ≈ $(round(std(lf_dist), digits=4)))")

# Plot histogram
histogram(lf_samples, bins=30, alpha=0.6, normalize=:pdf,
          xlabel="Load Factor", ylabel="Density",
          title="Triangular Distribution Test (mode=$(round(lf_mode, digits=2)))",
          legend=false)
savefig("_output/test_triangular.png")

println("\n✓ Triangular distribution test complete. Plot saved to _output/test_triangular.png")

println("\n"*"="^80)
println("TEST 3: Full gen_rand_vars() Integration Test")
println("="^80)

# Run gen_rand_vars on first project to verify everything works together
test_project = pjs[1]
electricity_price_mean = mean([52.2, 95.8])

println("\nTesting gen_rand_vars() with project: $(test_project.name)")
println("  Scale: $(test_project.scale)")
println("  Type: $(test_project.type)")
println("  Capacity: $(test_project.plant_capacity) MWe")

rand_vars = gen_rand_vars(
    "rothwell", n_test, wacc, electricity_price_mean, test_project;
    construction_time_range=ct_range
)

println("\nGenerated random variables:")
println("  WACC:")
println("    Range: [$(round(minimum(rand_vars.wacc), digits=4)), $(round(maximum(rand_vars.wacc), digits=4))]")
println("    Mean: $(round(mean(rand_vars.wacc), digits=4))")
println("  Construction Time:")
println("    Range: [$(minimum(rand_vars.construction_time)), $(maximum(rand_vars.construction_time))] years")
println("    Mean: $(round(mean(rand_vars.construction_time), digits=2)) years")
println("  Load Factor:")
println("    Range: [$(round(minimum(rand_vars.loadfactor), digits=4)), $(round(maximum(rand_vars.loadfactor), digits=4))]")
println("    Mean: $(round(mean(rand_vars.loadfactor), digits=4))")
println("  Investment:")
println("    Range: [$(round(minimum(rand_vars.investment)/1e6, digits=2)), $(round(maximum(rand_vars.investment)/1e6, digits=2))] M USD")
println("    Mean: $(round(mean(rand_vars.investment)/1e6, digits=2)) M USD")

# Verify correlation
actual_corr_full = cor(rand_vars.wacc, Float64.(rand_vars.construction_time))
println("\n  WACC × Construction Time correlation: $(round(actual_corr_full, digits=3))")
println("  Expected: ≈ 0.4")

if abs(actual_corr_full - 0.4) < 0.05
    println("  ✓ Correlation within ±0.05 of target")
else
    println("  ⚠ Correlation outside expected range")
end

println("\n✓ Integration test complete")

println("\n"*"="^80)
println("ALL TESTS COMPLETE")
println("="^80)
println("\nNext steps:")
println("  1. Review plots in _output/ directory")
println("  2. Verify correlation is within ±0.05 of 0.4")
println("  3. Check that distributions are within expected bounds")
println("  4. If tests pass, implement Shapley sensitivity")
println("="^80)
