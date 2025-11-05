"""
Simple script to plot OCC vs Year for Large reactors only
Shows Western (blue) vs Chinese/Asian (orange) with linear trend lines

Usage:
  julia plot_large_occ_vs_year.jl

Requirements:
- Run convert_to_csv.jl first to generate reactor_data.csv with year column
"""

using CairoMakie
using CSV
using DataFrames
using Statistics
using GLM  # For linear regression

# Set paths
inputpath = "_input"
outputpath = "_output"

println("="^60)
println("OCC vs Year Plot (Large Reactors Only)")
println("="^60)
println()

# Check if CSV file exists
reactor_csv_path = "$inputpath/reactor_data.csv"
if !isfile(reactor_csv_path)
    println("❌ ERROR: reactor_data.csv not found at $reactor_csv_path")
    println("Please run: julia convert_to_csv.jl first")
    exit(1)
end

# Read CSV
println("Reading reactor data...")
df = CSV.File(reactor_csv_path) |> DataFrame
println("✓ Loaded $(nrow(df)) total reactors")

# Check for year column
if !(:year in propertynames(df))
    println("❌ ERROR: year column not found in reactor_data.csv")
    println("Please regenerate the CSV by running: julia convert_to_csv.jl")
    exit(1)
end

# Filter to Large reactors only
df_large = filter(row -> row.scale == "Large", df)
println("✓ Filtered to $(nrow(df_large)) Large reactors")

if nrow(df_large) == 0
    println("❌ ERROR: No Large reactors found in dataset")
    exit(1)
end

# Calculate OCC (USD/kW) from investment (USD/MW)
df_large.occ = df_large.investment ./ 1000

# Group by region: Western vs Asian
# Define region groups
western_regions = ["Western / Developed", "Eastern Europe", "South America"]
asian_regions = ["Emerging Asia"]  # This includes China, India, Pakistan, Russia

df_large.group = map(df_large.region) do reg
    if reg in western_regions
        return "Western"
    elseif reg in asian_regions
        return "Asian"
    else
        return "Other"
    end
end

# Filter out "Other" if any
df_large = filter(row -> row.group != "Other", df_large)
println()
println("--- Regional Grouping ---")
println("Western reactors: ", count(df_large.group .== "Western"))
println("Asian reactors: ", count(df_large.group .== "Asian"))
println()

# Create the plot
fig = Figure(size=(900, 600))
ax = Axis(fig[1, 1],
          xlabel="Grid Connection Year",
          ylabel="OCC (USD2020/kW)",
          title="Overnight Capital Cost vs Year (Large Reactors)")

# Plot Western reactors (blue)
df_western = filter(row -> row.group == "Western", df_large)
if nrow(df_western) > 0
    scatter!(ax, df_western.year, df_western.occ,
            markersize=12, color=:blue, strokewidth=1, strokecolor=:black,
            label="Western")

    # Add linear trend line for Western
    if nrow(df_western) > 1
        # Fit linear model: OCC ~ Year
        lm_western = lm(@formula(occ ~ year), df_western)

        # Get coefficients
        coef_western = coef(lm_western)
        intercept_w = coef_western[1]
        slope_w = coef_western[2]

        # Generate trend line
        year_range_w = minimum(df_western.year):maximum(df_western.year)
        occ_pred_w = intercept_w .+ slope_w .* year_range_w

        lines!(ax, year_range_w, occ_pred_w,
              color=:blue, linewidth=2, linestyle=:dash,
              label="Western Trend")

        println("Western trend: OCC = $(round(intercept_w, digits=0)) + $(round(slope_w, digits=1)) × Year")
    end
end

# Plot Asian reactors (orange)
df_asian = filter(row -> row.group == "Asian", df_large)
if nrow(df_asian) > 0
    scatter!(ax, df_asian.year, df_asian.occ,
            markersize=12, color=:orange, strokewidth=1, strokecolor=:black,
            label="Asian")

    # Add linear trend line for Asian
    if nrow(df_asian) > 1
        # Fit linear model: OCC ~ Year
        lm_asian = lm(@formula(occ ~ year), df_asian)

        # Get coefficients
        coef_asian = coef(lm_asian)
        intercept_a = coef_asian[1]
        slope_a = coef_asian[2]

        # Generate trend line
        year_range_a = minimum(df_asian.year):maximum(df_asian.year)
        occ_pred_a = intercept_a .+ slope_a .* year_range_a

        lines!(ax, year_range_a, occ_pred_a,
              color=:orange, linewidth=2, linestyle=:dash,
              label="Asian Trend")

        println("Asian trend: OCC = $(round(intercept_a, digits=0)) + $(round(slope_a, digits=1)) × Year")
    end
end

# Add legend
axislegend(ax, position=:rt)

# Add grid
ax.xgridvisible = true
ax.ygridvisible = true

# Create output directory if needed
if !isdir(outputpath)
    mkdir(outputpath)
end

# Save plot
output_file = "$outputpath/fig-occ_vs_year-large_reactors.pdf"
save(output_file, fig)

println()
println("="^60)
println("✅ Plot saved to: $output_file")
println("="^60)
