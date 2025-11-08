#!/usr/bin/env julia

"""
Test script to verify reactor-type-specific capacity factor implementation
"""

println("="^70)
println("Testing Reactor-Type-Specific Capacity Factor Ranges")
println("="^70)

# Load conversion functions
include("convert_to_csv.jl")

# Test the capacity factor range function for all reactor types
test_types = ["BWR", "PWR", "HTR", "SFR", "LFR", "MSR", "MSFR", "MR", "HTGR", "HTR/GFR"]

println("\n1. Testing get_capacity_factor_range() function:")
println("-" * "70")
println("Reactor Type | CF Min | CF Max | Range Width | Midpoint")
println("-" * "70")

for rtype in test_types
    try
        cf_min, cf_max = get_capacity_factor_range(rtype)
        range_width = cf_max - cf_min
        midpoint = (cf_min + cf_max) / 2
        @printf("%-12s | %5.2f  | %5.2f  | %11.2f | %8.2f\n",
                rtype, cf_min, cf_max, range_width, midpoint)
    catch e
        println("$rtype: ERROR - $e")
    end
end

# Regenerate CSV with new capacity factor logic
println("\n2. Regenerating reactor_data.csv with type-specific capacity factors...")
println("-" * "70")

excel_file = "_input/reactor_data_raw.xlsx"
if isfile(excel_file)
    df = extract_reactor_data(excel_file)

    # Show capacity factors by reactor type
    println("\n3. Capacity factor ranges in generated CSV:")
    println("-" * "70")
    println("Reactor Name                | Type | CF Lower | CF Upper")
    println("-" * "70")

    for i in 1:min(20, nrow(df))  # Show first 20 reactors
        @printf("%-27s | %-4s | %8.3f | %8.3f\n",
                df.name[i], df.type[i], df.loadfactor_lower[i], df.loadfactor_upper[i])
    end

    # Summary statistics by type
    println("\n4. Summary by reactor type:")
    println("-" * "70")

    using Statistics
    for rtype in unique(df.type)
        mask = df.type .== rtype
        count = sum(mask)
        cf_min = unique(df.loadfactor_lower[mask])[1]
        cf_max = unique(df.loadfactor_upper[mask])[1]

        @printf("%-10s: %2d reactors → CF range [%.2f, %.2f]\n",
                rtype, count, cf_min, cf_max)
    end

    # Verify specific test cases from user requirements
    println("\n5. Validation of user requirements:")
    println("-" * "70")

    # Test BWRX-300 (BWR)
    bwr_mask = df.name .== "BWRX-300"
    if any(bwr_mask)
        bwr_cf = (df.loadfactor_lower[bwr_mask][1], df.loadfactor_upper[bwr_mask][1])
        expected_bwr = (0.75, 0.95)
        match_bwr = bwr_cf == expected_bwr
        println("✓ BWRX-300 (BWR): $(bwr_cf) $(match_bwr ? "✓ CORRECT" : "✗ MISMATCH")")
        println("  Expected: $(expected_bwr)")
    end

    # Test fast reactors (SFR)
    sfr_mask = df.type .== "SFR"
    if any(sfr_mask)
        sfr_cf = (df.loadfactor_lower[sfr_mask][1], df.loadfactor_upper[sfr_mask][1])
        expected_sfr = (0.55, 0.85)
        match_sfr = sfr_cf == expected_sfr
        println("✓ Fast reactors (SFR): $(sfr_cf) $(match_sfr ? "✓ CORRECT" : "✗ MISMATCH")")
        println("  Expected: $(expected_sfr)")
    end

    # Test PWR reactors
    pwr_mask = df.type .== "PWR"
    if any(pwr_mask)
        pwr_cf = (df.loadfactor_lower[pwr_mask][1], df.loadfactor_upper[pwr_mask][1])
        expected_pwr = (0.65, 0.95)
        match_pwr = pwr_cf == expected_pwr
        println("✓ PWR reactors: $(pwr_cf) $(match_pwr ? "✓ CORRECT" : "✗ MISMATCH")")
        println("  Expected: $(expected_pwr)")
    end

    println("\n" * "="^70)
    println("✓ Test completed successfully!")
    println("="^70)
else
    println("ERROR: Excel file not found: $excel_file")
    println("Please ensure reactor_data_raw.xlsx exists in _input/ directory")
end
