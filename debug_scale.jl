using XLSX
using DataFrames

println("Checking for Scale column in Excel file...")
println("="^60)

excel_file = "_input/reactor_data_raw.xlsx"
xf = XLSX.readxlsx(excel_file)
sheet = xf["Datasheet"]
df = DataFrame(XLSX.gettable(sheet; first_row=3))

println("\nAll column names:")
for (i, col) in enumerate(names(df))
    println("$i. \"$col\"")
end

println("\n" * "="^60)
println("Looking for Scale/scale columns specifically:")
println("="^60)

for col in names(df)
    if occursin(r"scale"i, col) || occursin(r"large.*medium.*micro"i, col)
        println("✓ Found: \"$col\"")
        println("  First 5 values: ", df[1:min(5, nrow(df)), col])
    end
end

println("\n" * "="^60)
println("Now checking the generated CSV file...")
println("="^60)

if isfile("_input/reactor_data.csv")
    using CSV
    csv_df = CSV.File("_input/reactor_data.csv"; delim=';') |> DataFrame
    println("\nCSV column names:")
    for (i, col) in enumerate(names(csv_df))
        println("$i. $col")
    end

    if "scale" in names(csv_df)
        println("\n✓ Scale column found in CSV!")
        println("First 5 values: ", csv_df[1:min(5, nrow(csv_df)), "scale"])
    else
        println("\n✗ Scale column NOT found in CSV")
    end
else
    println("\n✗ reactor_data.csv doesn't exist yet. Run julia convert_to_csv.jl first")
end
