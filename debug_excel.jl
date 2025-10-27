using XLSX
using DataFrames

println("="^80)
println("DEBUG: Excel File Column Inspector")
println("="^80)

excel_file = "_input/reactor_data_raw.xlsx"

# Read the Excel file
xf = XLSX.readxlsx(excel_file)
sheet = xf["Datasheet"]

println("\nReading with first_row=3 (skiprows=2)...")
df = DataFrame(XLSX.gettable(sheet; first_row=3))

println("\nNumber of rows: ", nrow(df))
println("Number of columns: ", ncol(df))

println("\n" * "="^80)
println("ALL COLUMN NAMES:")
println("="^80)
for (i, col) in enumerate(names(df))
    println("$i. \"$col\"")
end

println("\n" * "="^80)
println("FIRST 3 DATA ROWS (selected columns):")
println("="^80)

# Show key columns
key_cols = []
for col in names(df)
    if occursin(r"Project|Type|Capacity|Construction|Lifetime|OPEX|Fuel"i, col)
        push!(key_cols, col)
    end
end

println("\nKey columns found:")
for col in key_cols
    println("  - $col")
end

println("\n" * "="^80)
println("SAMPLE DATA FROM KEY COLUMNS:")
println("="^80)

# Print specific columns
if "Project" in names(df)
    println("\nProject (first 5):")
    println(df[1:min(5,nrow(df)), "Project"])
end

if "planned Construction time (y)" in names(df)
    println("\nPlanned Construction time (first 5):")
    println(df[1:min(5,nrow(df)), "planned Construction time (y)"])
end

if "Actual / total build time (y)" in names(df)
    println("\nActual/total build time (first 5):")
    println(df[1:min(5,nrow(df)), "Actual / total build time (y)"])
end

if "Lifetime (y)" in names(df)
    println("\nLifetime (first 5):")
    println(df[1:min(5,nrow(df)), "Lifetime (y)"])
end

if "Type" in names(df)
    println("\nType (first 5):")
    println(df[1:min(5,nrow(df)), "Type"])
end

println("\n" * "="^80)
println("CHECKING FOR DUPLICATED HEADER ROWS:")
println("="^80)

if "Type" in names(df)
    type_col = df[!, "Type"]
    for (i, val) in enumerate(type_col[1:min(10, length(type_col))])
        println("Row $i Type value: \"$val\" (type: $(typeof(val)))")
    end
end
