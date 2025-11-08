using Pkg
Pkg.activate(pwd())
using Statistics
using Logging

# Suppress @info logging to reduce terminal spam
global_logger(SimpleLogger(stdout, Logging.Warn))

# Define paths required by data.jl
inputpath = "_input"
outputpath = "_output"

include("functions.jl")
include("data.jl")

# Test parameters
n_test = 500
n_outer = 50  # Define early so we can use it in print
n_inner = 100
wacc = [0.04, 0.10]
ct_range = [3, 7]
scaling = [0.20, 0.75]

# Open results file
io = open("diagnostic_results.txt", "w")

println(io, "="^80)
println(io, "SHAPLEY CONDITIONAL VARIANCE DIAGNOSTICS")
println(io, "="^80)

# Select test project
test_project = pjs[1]
electricity_price_mean = mean([52.2, 95.8])

println(io, "\nTest Project: $(test_project.name)")
println(io, "  Capacity: $(test_project.plant_capacity) MWe")
println(io, "  Parameters: n_outer=$n_outer, n_inner=$n_inner, n_test=$n_test")

# Seed and generate A/B samples
using Random
Random.seed!(54321)
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

println(io, "\n" * "="^80)
println(io, "TOTAL VARIANCE (from A/B samples)")
println(io, "="^80)
println(io, "  NPV:  $(round(total_var_npv, sigdigits=6))")
println(io, "  LCOE: $(round(total_var_lcoe, sigdigits=6))")

# Generate outer base (n_outer and n_inner defined at top of file)
Random.seed!(12345)
outer_base = []
for k in 1:n_outer
    sample = gen_rand_vars("rothwell", 1, wacc, electricity_price_mean, test_project;
                          construction_time_range=ct_range)
    push!(outer_base, sample)
end

println(io, "\n" * "="^80)
println(io, "TEST 1: V(∅) - Empty Coalition")
println(io, "="^80)
println(io, "Expected: < 1% of total variance")

coalition_id_empty = hash(Symbol[])
V_empty_npv, V_empty_lcoe = estimate_conditional_variance_V_S(
    test_project, Symbol[], outer_base, n_inner,
    "rothwell", wacc, electricity_price_mean,
    ct_range, coalition_id_empty
)

ratio_npv_empty = V_empty_npv/total_var_npv
ratio_lcoe_empty = V_empty_lcoe/total_var_lcoe

println(io, "\nNPV:")
println(io, "  V(∅):        $(round(V_empty_npv, sigdigits=6))")
println(io, "  Total var:   $(round(total_var_npv, sigdigits=6))")
println(io, "  Ratio:       $(round(ratio_npv_empty, digits=4)) (MUST be < 0.01)")
println(io, "  Status:      $(ratio_npv_empty < 0.01 ? "✓ PASS" : "❌ FAIL")")

println(io, "\nLCOE:")
println(io, "  V(∅):        $(round(V_empty_lcoe, sigdigits=6))")
println(io, "  Total var:   $(round(total_var_lcoe, sigdigits=6))")
println(io, "  Ratio:       $(round(ratio_lcoe_empty, digits=4)) (MUST be < 0.01)")
println(io, "  Status:      $(ratio_lcoe_empty < 0.01 ? "✓ PASS" : "❌ FAIL")")

println(io, "\n" * "="^80)
println(io, "TEST 2: V(all) - Full Coalition")
println(io, "="^80)
println(io, "Expected: 90-110% of total variance")

coalition_id_full = hash([:wacc, :construction_time, :loadfactor, :investment])
V_all_npv, V_all_lcoe = estimate_conditional_variance_V_S(
    test_project, [:wacc, :construction_time, :loadfactor, :investment],
    outer_base, n_inner, "rothwell", wacc, electricity_price_mean,
    ct_range, coalition_id_full
)

ratio_npv_all = V_all_npv/total_var_npv
ratio_lcoe_all = V_all_lcoe/total_var_lcoe

println(io, "\nNPV:")
println(io, "  V(all):      $(round(V_all_npv, sigdigits=6))")
println(io, "  Total var:   $(round(total_var_npv, sigdigits=6))")
println(io, "  Ratio:       $(round(ratio_npv_all, digits=4)) (MUST be 0.9-1.1)")
println(io, "  Status:      $(0.9 <= ratio_npv_all <= 1.1 ? "✓ PASS" : "❌ FAIL")")

println(io, "\nLCOE:")
println(io, "  V(all):      $(round(V_all_lcoe, sigdigits=6))")
println(io, "  Total var:   $(round(total_var_lcoe, sigdigits=6))")
println(io, "  Ratio:       $(round(ratio_lcoe_all, digits=4)) (MUST be 0.9-1.1)")
println(io, "  Status:      $(0.9 <= ratio_lcoe_all <= 1.1 ? "✓ PASS" : "❌ FAIL")")

println(io, "\n" * "="^80)
println(io, "SUMMARY")
println(io, "="^80)

test1_pass = ratio_npv_empty < 0.01 && ratio_lcoe_empty < 0.01
test2_pass = (0.9 <= ratio_npv_all <= 1.1) && (0.9 <= ratio_lcoe_all <= 1.1)

println(io, "\nTest 1 (V(∅) ≈ 0):     $(test1_pass ? "✓ PASS" : "❌ FAIL")")
println(io, "Test 2 (V(all) ≈ var): $(test2_pass ? "✓ PASS" : "❌ FAIL")")

if !test1_pass
    println(io, "\n⚠️  CRITICAL: V(∅) captures $(round(100*ratio_lcoe_empty, digits=1))% of variance!")
    println(io, "    This means the empty coalition isn't varying all parameters.")
    println(io, "    Bug likely in: generate_combined_sample_from_outer_base()")
end

close(io)

# Print to terminal
println("\n" * "="^80)
println("DIAGNOSTIC TEST COMPLETE")
println("="^80)
println("Results saved to: diagnostic_results.txt")
println("\nQuick Summary:")
println("  Test 1 (V(∅) ≈ 0):     $(test1_pass ? "✓ PASS" : "❌ FAIL - V(∅) ratio = $(round(ratio_lcoe_empty, digits=3))")")
println("  Test 2 (V(all) ≈ var): $(test2_pass ? "✓ PASS" : "❌ FAIL - V(all) ratio = $(round(ratio_lcoe_all, digits=3))")")
println("\nSee diagnostic_results.txt for full details.")
println("="^80)
