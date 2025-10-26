# Reactor Data Conversion Instructions

This document explains how to convert the Excel reactor data to CSV format using Julia.

## Prerequisites

You need to have Julia installed with the required packages.

## Installing Required Packages

### Option 1: Using the Julia REPL (Recommended)

1. Open Julia in your project directory
2. Enter package mode by pressing `]`
3. Run these commands:

```julia
activate .
add XLSX
add DataFrames
add CSV
instantiate
```

Then press backspace to exit package mode.

### Option 2: Command Line (Quick)

Run this command from your terminal in the project directory:

```bash
julia --project=. -e "using Pkg; Pkg.add(\"XLSX\"); Pkg.add(\"DataFrames\"); Pkg.add(\"CSV\"); Pkg.instantiate()"
```

### Option 3: From Julia REPL (Alternative)

```julia
using Pkg
Pkg.activate(".")
Pkg.add("XLSX")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.instantiate()
```

## Running the Conversion Script

Once packages are installed, run:

```bash
julia convert_to_csv.jl
```

Or from within Julia:

```julia
include("convert_to_csv.jl")
```

## Output

The script will:
- Read `_input/reactor_data_raw.xlsx`
- Process and standardize reactor types
- Convert to the same CSV format as `project_data.csv`
- Save output to `_input/reactor_data.csv`

## CSV Structure

The output CSV will have these columns (semicolon-delimited):
1. `name` - Reactor project name
2. `type` - Standardized reactor type (PWR, BWR, HTR, SFR, etc.)
3. `investment` - Total investment in USD
4. `plant_capacity` - Capacity in MWe
5. `learning_factor` - Learning rate (default: 0.1)
6. `construction_time` - Construction time in years
7. `operating_time` - Operating lifetime in years
8. `loadfactor_lower` - Lower bound of capacity factor (decimal)
9. `loadfactor_upper` - Upper bound of capacity factor (decimal)
10. `operating_cost_fix` - Fixed operating costs (USD/year)
11. `operating_cost_variable` - Variable operating costs (USD/MWh)
12. `operating_cost_fuel` - Fuel costs (USD/MWh)
13. `reference_pj_investment` - Reference project investment
14. `reference_pj_capacity` - Reference project capacity

## Troubleshooting

### Error: "Package XLSX not found"
- You need to install the packages first (see above)

### Error: "UndefVarError: `Pkg` not defined"
- You need to load Pkg first: `using Pkg`

### Error: "Cannot open Excel file"
- Make sure `_input/reactor_data_raw.xlsx` exists
- Check that you're running from the project root directory
