using XLSX
using DataFrames
using CSV
using Statistics

"""
    standardize_reactor_type(reactor_type)

Standardize reactor type names to basic technology categories

# Arguments
- `reactor_type::String`: Original reactor type from data

# Returns
- `String`: Standardized reactor type
"""
function standardize_reactor_type(reactor_type)
    # Convert to uppercase and strip whitespace
    rtype = uppercase(strip(string(reactor_type)))

    # PWR variants
    pwr_types = [
        "PWR", "EPR", "EPR (PWR)", "EPR-1750",
        "AP1000", "APR1400",
        "VVER1200 (PWR)", "VVER", "VVER1200",
        "IPWR"  # Integral PWR
    ]

    # BWR variants
    bwr_types = ["BWR"]

    # HTR variants
    htr_types = ["HTR", "HTR/GFR", "HTGR"]

    # SFR variants
    sfr_types = ["SFR", "BN-800 (SFR)", "BN-800", "BN800"]

    # MSR variants
    msr_types = ["MSR", "MSFR"]

    # LFR variants
    lfr_types = ["LFR"]

    # MR variants
    mr_types = ["MR"]

    # Check which category it belongs to
    if any(occursin(pwr, rtype) for pwr in pwr_types)
        return "PWR"
    elseif any(occursin(bwr, rtype) for bwr in bwr_types)
        return "BWR"
    elseif any(occursin(htr, rtype) for htr in htr_types)
        return "HTR"
    elseif any(occursin(sfr, rtype) for sfr in sfr_types)
        return "SFR"
    elseif any(occursin(msr, rtype) for msr in msr_types)
        return "MSR"
    elseif any(occursin(lfr, rtype) for lfr in lfr_types)
        return "LFR"
    elseif any(occursin(mr, rtype) for mr in mr_types)
        return "MR"
    else
        println("⚠ Warning: Unknown reactor type '$reactor_type', keeping as-is")
        return reactor_type
    end
end


"""
    assign_tech_category(size_category, reactor_type)

Assign combined technology category (Size-Type)

# Arguments
- `size_category::String`: Size category (SMR, Large, etc.)
- `reactor_type::String`: Reactor type (PWR, BWR, etc.)

# Returns
- `String`: Combined category like 'SMR-PWR', 'Large-BWR', etc.
"""
function assign_tech_category(size_category, reactor_type)
    size = ismissing(size_category) ? "Unknown" : size_category
    rtype = ismissing(reactor_type) ? "Unknown" : reactor_type
    return "$size-$rtype"
end


"""
    clean_numeric_string(s)

Clean numeric strings by removing spaces, dashes, and 'x' markers
"""
function clean_numeric_string(s)
    if ismissing(s)
        return missing
    end

    str = string(s)
    str = replace(str, " " => "")
    str = replace(str, "–" => "")
    str = replace(str, "-" => "")
    str = replace(str, "x" => "")

    if isempty(str)
        return missing
    end

    try
        return parse(Float64, str)
    catch
        return missing
    end
end


"""
    extract_reactor_data(excel_file::String; output_csv::String="_input/reactor_data.csv")

Extract reactor data from Excel for LCOE Monte Carlo simulation

# Arguments
- `excel_file::String`: Path to your Excel file
- `output_csv::String`: Output CSV filename (default: "_input/reactor_data.csv")

# Returns
- `DataFrame`: Processed reactor data
"""
function extract_reactor_data(excel_file::String; output_csv::String="_input/reactor_data.csv")
    # Read the specific sheet
    xf = XLSX.readxlsx(excel_file)
    sheet = xf["Datasheet"]

    # Read data starting from row 3 (skip first 2 rows)
    df = DataFrame(XLSX.gettable(sheet; first_row=3))

    println("Available columns in Excel:")
    for (i, col) in enumerate(names(df))
        println("$i. '$col'")
    end
    println()

    # Define column mapping for LCOE simulation
    column_mapping = Dict(
        # Identifiers
        "Country" => :country,
        "Project" => :name,
        "Type" => :reactor_type,
        "FOAK/NOAK*" => :foak_noak,
        "Large medium micro" => :size_category,

        # Core economic parameters
        "Capacity (net MWe)" => :capacity_mwe,
        "OCC (USD2020/kW)" => :occ_usd_per_kw,

        # Construction time parameters
        "planned Construction time (y)" => :construction_time_planned,
        "Actual / total build time (y)" => :construction_time_actual,

        # Operational parameters
        "Lifetime (y)" => :lifetime_years,
        "Capacity factor (%)" => :capacity_factor_pct,

        # Cost parameters
        "OPEX fixed (USD2020/MW-yr)" => :opex_fixed_usd_per_mw_yr,
        "OPEX variable (USD2020/MWh)" => :opex_variable_usd_per_mwh,
        "Fuel (USD2020/MWh)" => :fuel_usd_per_mwh,
        "Waste (USD2020/MWh)" => :waste_usd_per_mwh,
        "Decom (% of CAPEX)" => :decom_pct_of_capex,

        # Optional but useful
        "WACC (real %)" => :wacc_pct,
        "Year" => :data_year
    )

    # Check which columns exist
    available_cols = Dict{String, Symbol}()
    missing_cols = String[]

    for (old_name, new_name) in column_mapping
        if old_name in names(df)
            available_cols[old_name] = new_name
        else
            push!(missing_cols, old_name)
        end
    end

    if !isempty(missing_cols)
        println("Warning: These columns not found: $missing_cols\n")
    end

    if isempty(available_cols)
        println("ERROR: No matching columns found!")
        return nothing
    end

    # Select and rename columns
    df_output = select(df, keys(available_cols)...)
    rename!(df_output, available_cols)

    # Clean numeric columns
    numeric_cols = [
        :capacity_mwe, :occ_usd_per_kw,
        :construction_time_planned, :construction_time_actual,
        :lifetime_years, :capacity_factor_pct,
        :opex_fixed_usd_per_mw_yr, :opex_variable_usd_per_mwh,
        :fuel_usd_per_mwh, :waste_usd_per_mwh,
        :decom_pct_of_capex, :wacc_pct
    ]

    for col in numeric_cols
        if col in names(df_output)
            df_output[!, col] = clean_numeric_string.(df_output[!, col])
        end
    end

    # === Standardize reactor types ===
    println("--- Reactor Type Standardization ---")
    if :reactor_type in names(df_output)
        println("Original reactor types:")
        println(combine(groupby(df_output, :reactor_type), nrow => :count))

        df_output[!, :reactor_type] = standardize_reactor_type.(df_output[!, :reactor_type])

        println("\nStandardized reactor types:")
        println(combine(groupby(df_output, :reactor_type), nrow => :count))
        println()
    end

    # === Assign combined tech category ===
    if :size_category in names(df_output) && :reactor_type in names(df_output)
        df_output[!, :tech_category] = assign_tech_category.(df_output[!, :size_category], df_output[!, :reactor_type])
        println("--- Technology Categories Created ---")
        println(combine(groupby(df_output, :tech_category), nrow => :count))
        println()
    end

    # Data validation and cleaning
    println("--- Data Cleaning ---")

    # Calculate absolute OCC (total capital cost)
    if :occ_usd_per_kw in names(df_output) && :capacity_mwe in names(df_output)
        df_output[!, :occ_total_million_usd] = (
            df_output[!, :occ_usd_per_kw] .* df_output[!, :capacity_mwe] .* 1000 ./ 1e6
        )
    end

    # Use actual construction time if available, otherwise planned
    if :construction_time_actual in names(df_output)
        if :construction_time_planned in names(df_output)
            df_output[!, :construction_time] = coalesce.(df_output[!, :construction_time_actual], df_output[!, :construction_time_planned])
        else
            df_output[!, :construction_time] = df_output[!, :construction_time_actual]
        end
    end

    # Remove rows with missing critical data
    critical_cols = [:name, :capacity_mwe, :occ_usd_per_kw]
    existing_critical = [col for col in critical_cols if col in names(df_output)]

    if !isempty(existing_critical)
        df_output = dropmissing(df_output, existing_critical)
    end

    # Sort by tech category for easier analysis
    if :tech_category in names(df_output)
        sort!(df_output, :tech_category)
    end

    # Save as semicolon-delimited CSV
    CSV.write(output_csv, df_output; delim=';')

    println("\n✓ Successfully extracted $(nrow(df_output)) reactors")
    println("✓ Saved to: $output_csv")

    # Summary by tech category
    if :tech_category in names(df_output)
        println("\n--- Reactors by Technology Category ---")
        println(combine(groupby(df_output, :tech_category), nrow => :count))

        println("\n--- Cost Range by Category (USD/kW) ---")
        cost_summary = combine(groupby(df_output, :tech_category)) do sdf
            (
                count = nrow(sdf),
                min = minimum(skipmissing(sdf.occ_usd_per_kw)),
                median = median(skipmissing(sdf.occ_usd_per_kw)),
                max = maximum(skipmissing(sdf.occ_usd_per_kw))
            )
        end
        println(cost_summary)
    end

    # Summary by size category
    if :size_category in names(df_output)
        println("\n--- Reactors by Size Category ---")
        println(combine(groupby(df_output, :size_category), nrow => :count))
    end

    # Summary by reactor type
    if :reactor_type in names(df_output)
        println("\n--- Reactors by Type ---")
        println(combine(groupby(df_output, :reactor_type), nrow => :count))
    end

    println("\n--- Columns in Output ---")
    println(names(df_output))

    println("\n--- Data Completeness ---")
    for col in names(df_output)
        missing_count = count(ismissing, df_output[!, col])
        if missing_count > 0
            println("$col: $missing_count missing")
        end
    end

    return df_output
end


# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    excel_file = "_input/reactor_data_raw.xlsx"
    df = extract_reactor_data(excel_file)

    if !isnothing(df)
        println("\n--- First Few Rows ---")
        display_cols = [:name, :reactor_type, :size_category, :tech_category, :occ_usd_per_kw]
        existing_display_cols = [col for col in display_cols if col in names(df)]
        println(first(df[!, existing_display_cols], 10))

        # Check for FOAK vs NOAK distribution
        if :foak_noak in names(df) && :tech_category in names(df)
            println("\n--- FOAK/NOAK Distribution by Tech Category ---")
            println(combine(groupby(df, [:tech_category, :foak_noak]), nrow => :count))
        end
    end
end
