#!/bin/bash
# Test script to verify convert_to_csv.jl works

echo "Testing reactor data conversion..."
echo ""

# Create a simple Julia script that calls the function
cat > /tmp/test_conversion.jl <<'EOF'
include("convert_to_csv.jl")
df = extract_reactor_data("_input/reactor_data_raw.xlsx")
println("\n=== TEST COMPLETE ===")
println("Dataframe has $(nrow(df)) rows and $(ncol(df)) columns")
println("Column names: ", names(df))
EOF

# Try to run it
if command -v julia &> /dev/null; then
    julia /tmp/test_conversion.jl
    echo ""
    echo "CSV file line count:"
    wc -l _input/reactor_data.csv
else
    echo "Julia not available in this environment"
    echo "Please run this on your local machine:"
    echo "  julia convert_to_csv.jl"
fi
