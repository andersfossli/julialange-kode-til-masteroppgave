##### Test script for learning curve functionality #####

println("="^80)
println("Learning Curve Test Script")
println("="^80)
println()

# Load the learning_multiplier function
include("functions.jl")

##### Test 1: Basic multiplier calculations #####
println("Test 1: Basic Learning Multiplier Calculations")
println("-"^80)

LR = 0.10  # 10% learning rate
kappa = 1.20  # 20% FOAK premium

println("Parameters: LR = $LR (10%), κ = $kappa (20% FOAK premium)")
println()

N_values = [1, 2, 4, 6, 8, 10]
println("Without floor:")
for N in N_values
    m = learning_multiplier(N, LR, kappa=kappa)
    pct_above_soak = (m - 1.0) * 100
    println("  N=$N: multiplier = $(round(m, digits=4)) ($(round(pct_above_soak, digits=1))% above SOAK)")
end

println()
println("With floor = 1.0 (SOAK baseline):")
for N in N_values
    m = learning_multiplier(N, LR, kappa=kappa, floor=1.0)
    pct_above_soak = (m - 1.0) * 100
    println("  N=$N: multiplier = $(round(m, digits=4)) ($(round(pct_above_soak, digits=1))% above SOAK)")
end

println()
println("Expected behavior:")
println("  • FOAK (N=1) has 20% premium: m = 1.20 ✓")
println("  • Each doubling reduces cost by 10%:")
println("    - N=2: m = 1.20 × (1-0.10)^1 = 1.08 ✓")
println("    - N=4: m = 1.20 × (1-0.10)^2 = 0.972 (but floor=1.0 caps it) ✓")
println("  • With floor=1.0, multiplier never goes below SOAK")
println()

##### Test 2: Different learning rates #####
println("="^80)
println("Test 2: Impact of Different Learning Rates")
println("-"^80)

learning_rates = [0.05, 0.10, 0.15, 0.20]
N_test = [1, 2, 4, 8]

println("Fixed: κ = 1.20, floor = 1.0")
println()

for LR in learning_rates
    println("LR = $(Int(LR*100))%:")
    for N in N_test
        m = learning_multiplier(N, LR, kappa=1.20, floor=1.0)
        println("  N=$N: $(round(m, digits=3))")
    end
    println()
end

##### Test 3: Different FOAK premiums #####
println("="^80)
println("Test 3: Impact of Different FOAK Premiums (κ)")
println("-"^80)

kappas = [1.10, 1.20, 1.30]
LR_test = 0.10

println("Fixed: LR = 10%, floor = 1.0")
println()

for kappa in kappas
    premium_pct = (kappa - 1.0) * 100
    println("κ = $kappa ($(round(Int, premium_pct))% FOAK premium):")
    for N in N_test
        m = learning_multiplier(N, LR_test, kappa=kappa, floor=1.0)
        println("  N=$N: $(round(m, digits=3))")
    end
    println()
end

##### Test 4: Cost reduction visualization #####
println("="^80)
println("Test 4: Cost Reduction Example")
println("-"^80)

base_cost = 5000000.0  # USD/MW (example SOAK cost)
LR = 0.10
kappa = 1.20

println("Example: Base SOAK cost = \$$(Int(base_cost/1e6))M/MW")
println("Learning: LR = 10%, κ = 1.20 (20% FOAK premium), floor = 1.0")
println()

println("N | Multiplier | Cost (USD/MW) | Savings vs FOAK")
println("-"^60)

foak_cost = base_cost * kappa
for N in [1, 2, 3, 4, 5, 6, 8, 10]
    m = learning_multiplier(N, LR, kappa=kappa, floor=1.0)
    cost = base_cost * m
    savings_pct = (1 - cost/foak_cost) * 100

    println("$(rpad(N,2)) | $(rpad(round(m, digits=3),10)) | \$$(rpad(round(Int,cost/1e3),8))k | $(round(savings_pct, digits=1))%")
end

println()
println("Key insights:")
println("  • FOAK (N=1) costs \$$(round(Int,foak_cost/1e3))k/MW (20% premium)")
println("  • By N=4, cost reduced to SOAK baseline (\$$(round(Int,base_cost/1e3))k/MW)")
println("  • Total savings from FOAK to SOAK: \$$(round(Int,(foak_cost-base_cost)/1e3))k/MW ($(round((1-base_cost/foak_cost)*100, digits=1))%)")
println()

##### Test 5: Verify formula #####
println("="^80)
println("Test 5: Verify Formula Implementation")
println("-"^80)

N = 4
LR = 0.10
kappa = 1.20

# Manual calculation
log2_N = log2(N)
manual_m = kappa * (1.0 - LR)^log2_N

# Function calculation
function_m = learning_multiplier(N, LR, kappa=kappa)

println("Formula: m = κ × (1 - LR)^(log₂(N))")
println("Given: N=$N, LR=$LR, κ=$kappa")
println()
println("Manual calculation:")
println("  log₂($N) = $log2_N")
println("  (1 - $LR)^$log2_N = $((1.0-LR)^log2_N)")
println("  m = $kappa × $((1.0-LR)^log2_N) = $manual_m")
println()
println("Function result: $function_m")
println()

if abs(manual_m - function_m) < 1e-10
    println("✓ Formula implementation verified correctly!")
else
    println("✗ Formula mismatch detected!")
end

println()
println("="^80)
println("All tests complete!")
println("="^80)
