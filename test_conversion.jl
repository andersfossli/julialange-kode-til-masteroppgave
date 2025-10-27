println("Script started!")

# Test if packages are available
println("Testing package imports...")

try
    using XLSX
    println("✓ XLSX loaded successfully")
catch e
    println("✗ ERROR loading XLSX: ", e)
    println("\nYou need to install XLSX first:")
    println("  julia> using Pkg")
    println("  julia> Pkg.add(\"XLSX\")")
    exit(1)
end

try
    using DataFrames
    println("✓ DataFrames loaded successfully")
catch e
    println("✗ ERROR loading DataFrames: ", e)
    exit(1)
end

try
    using CSV
    println("✓ CSV loaded successfully")
catch e
    println("✗ ERROR loading CSV: ", e)
    exit(1)
end

using Statistics

println("\nAll packages loaded successfully!")
println("\nChecking for input file...")

excel_file = "_input/reactor_data_raw.xlsx"
if isfile(excel_file)
    println("✓ Found: $excel_file")
else
    println("✗ ERROR: Cannot find $excel_file")
    println("Current directory: ", pwd())
    println("Files in _input/:")
    if isdir("_input")
        for f in readdir("_input")
            println("  - $f")
        end
    else
        println("  _input/ directory not found!")
    end
    exit(1)
end

println("\nStarting conversion...")
println("="^60)

# Include the actual conversion code
include("convert_to_csv.jl")

println("\nCalling extract_reactor_data...")
df = extract_reactor_data(excel_file)

if !isnothing(df)
    println("\n✓ SUCCESS! Created reactor_data.csv with $(nrow(df)) reactors")
else
    println("\n✗ ERROR: Conversion failed")
end
