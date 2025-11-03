using XLSX
using DataFrames
using CSV
using Statistics

"""
    country_to_region(country)

Map country names to regional categories for analysis
"""
function country_to_region(country)
    if ismissing(country)
        return "Unknown"
    end

    country_str = uppercase(strip(string(country)))

    # East Asia
    east_asia = ["CHINA", "SOUTH KOREA", "KOREA", "JAPAN", "TAIWAN"]

    # Western
    western = ["USA", "UNITED STATES", "US", "FRANCE", "UK", "UNITED KINGDOM",
               "CANADA", "GERMANY", "SPAIN", "ITALY", "BELGIUM", "NETHERLANDS",
               "SWEDEN", "FINLAND", "SWITZERLAND"]

    # Eastern Europe / Russia
    eastern_europe = ["RUSSIA", "RUSSIAN FEDERATION", "UKRAINE", "CZECH REPUBLIC",
                      "SLOVAKIA", "BULGARIA", "ROMANIA", "HUNGARY", "POLAND"]

    # Middle East / South Asia
    middle_east_south_asia = ["INDIA", "PAKISTAN", "UAE", "SAUDI ARABIA",
                             "IRAN", "TURKEY", "EGYPT"]

    # South America
    south_america = ["BRAZIL", "ARGENTINA", "MEXICO"]

    # Check which region
    if any(occursin(c, country_str) for c in east_asia)
        return "East Asia"
    elseif any(occursin(c, country_str) for c in western)
        return "Western"
    elseif any(occursin(c, country_str) for c in eastern_europe)
        return "Eastern Europe"
    elseif any(occursin(c, country_str) for c in middle_east_south_asia)
        return "Middle East / South Asia"
    elseif any(occursin(c, country_str) for c in south_america)
        return "South America"
    else
        return "Other"
    end
end

"""
    standardize_reactor_type(reactor_type)

Standardize reactor type names to basic technology categories
"""
function standardize_reactor_type(reactor_type)
    rtype = uppercase(strip(string(reactor_type)))

    # PWR variants
    pwr_types = ["PWR", "EPR", "EPR (PWR)", "EPR-1750", "AP1000", "APR1400",
                 "VVER1200 (PWR)", "VVER", "VVER1200", "IPWR"]

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
        return string(reactor_type)
    end
end


"""
    get_reference_values(reactor_type)

Get reference project values based on reactor type (from existing data)
"""
function get_reference_values(reactor_type)
    if reactor_type == "PWR"
        return (investment=8600000, capacity=1117)
    elseif reactor_type == "BWR"
        return (investment=9722604, capacity=935)
    elseif reactor_type == "HTR"
        return (investment=7195125, capacity=330)
    elseif reactor_type == "SFR"
        return (investment=27747200, capacity=1250)
    elseif reactor_type == "MSR"
        return (investment=8600000, capacity=1000)  # Default, adjust as needed
    elseif reactor_type == "LFR"
        return (investment=27747200, capacity=1250)  # Similar to SFR
    else
        return (investment=8600000, capacity=1000)  # Default fallback
    end
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
    str = replace(str, "x" => "")

    # Handle ranges like "85-90" by taking the midpoint
    if occursin("-", str) && !startswith(str, "-")
        parts = split(str, "-")
        if length(parts) == 2
            try
                val1 = parse(Float64, parts[1])
                val2 = parse(Float64, parts[2])
                return (val1 + val2) / 2.0
            catch
                str = replace(str, "-" => "")
            end
        end
    else
        str = replace(str, "-" => "")
    end

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

Extract reactor data from Excel and format to match existing project_data.csv structure
"""
function extract_reactor_data(excel_file::String; output_csv::String="_input/reactor_data.csv")
    # Read the Excel file
    xf = XLSX.readxlsx(excel_file)
    sheet = xf["Datasheet"]
    df = DataFrame(XLSX.gettable(sheet; first_row=3))

    println("Available columns in Excel:")
    for (i, col) in enumerate(names(df))
        println("$i. '$col'")
    end
    println()

    # Define column mapping
    column_mapping = Dict(
        "Project" => :name,
        "Type" => :reactor_type_raw,
        "Country" => :country,  # New: for regional analysis
        "Capacity (net MWe)" => :capacity_mwe,
        "OCC (USD2020/kW)" => :occ_usd_per_kw,
        "planned Construction time (y)" => :construction_time_planned,
        "Actual / total build time (y)" => :construction_time_actual,
        "Lifetime (y)" => :lifetime_years,
        "Capacity factor (%)" => :capacity_factor_pct,
        "OPEX fixed (USD2020/MW-yr)" => :opex_fixed_usd_per_mw_yr,
        "OPEX variable (USD2020/MWh)" => :opex_variable_usd_per_mwh,
        "Fuel (USD2020/MWh)" => :fuel_usd_per_mwh,
        "Scale" => :scale,  # New column name (capital S)
        "scale" => :scale,  # New column name (lowercase s)
        "Large medium micro" => :scale,  # Old column name (backwards compatibility)
        "FOAK/NOAK*" => :foak_noak
    )

    # Check which columns exist
    available_cols = Dict{String, Symbol}()
    for (old_name, new_name) in column_mapping
        if old_name in names(df)
            available_cols[old_name] = new_name
        end
    end

    # Select and rename columns
    df_work = select(df, keys(available_cols)...)
    rename!(df_work, available_cols)

    # Convert all column names to symbols for consistency
    rename!(df_work, Symbol.(names(df_work)))

    # Filter out header rows that might have been included
    # Remove any row where the name or reactor_type_raw contains "Project" or "Type" (headers)
    if :name in propertynames(df_work)
        df_work = filter(row -> !ismissing(row.name) &&
                                !occursin(r"^Project$"i, string(row.name)), df_work)
    end
    if :reactor_type_raw in propertynames(df_work)
        df_work = filter(row -> !ismissing(row.reactor_type_raw) &&
                                !occursin(r"^Type$"i, string(row.reactor_type_raw)), df_work)
    end

    # Debug: Check which columns exist before cleaning
    println("\n--- Columns in df_work BEFORE cleaning ---")
    println("Column names: ", names(df_work))
    println("Property names (symbols): ", propertynames(df_work))
    if :construction_time_planned in propertynames(df_work)
        println("✓ construction_time_planned exists, sample: ", df_work[1:min(3, nrow(df_work)), :construction_time_planned])
    else
        println("✗ construction_time_planned NOT FOUND")
    end
    if :lifetime_years in propertynames(df_work)
        println("✓ lifetime_years exists, sample: ", df_work[1:min(3, nrow(df_work)), :lifetime_years])
    else
        println("✗ lifetime_years NOT FOUND")
    end

    # Clean numeric columns
    numeric_cols = [:capacity_mwe, :occ_usd_per_kw, :construction_time_planned,
                    :construction_time_actual, :lifetime_years, :capacity_factor_pct,
                    :opex_fixed_usd_per_mw_yr, :opex_variable_usd_per_mwh, :fuel_usd_per_mwh]

    println("\n--- Cleaning numeric columns ---")
    for col in numeric_cols
        if col in propertynames(df_work)
            println("Cleaning column: $col")
            before_sample = df_work[1:min(3, nrow(df_work)), col]
            println("  Before: $before_sample (types: $(typeof.(before_sample)))")
            df_work[!, col] = clean_numeric_string.(df_work[!, col])
            after_sample = df_work[1:min(3, nrow(df_work)), col]
            println("  After: $after_sample (types: $(typeof.(after_sample)))")
        end
    end

    # Remove rows with missing critical data
    initial_rows = nrow(df_work)
    df_work = dropmissing(df_work, [:name, :capacity_mwe, :occ_usd_per_kw])
    println("Rows after removing missing critical data: $(nrow(df_work)) (removed $(initial_rows - nrow(df_work)))")

    # Ensure numeric columns are actually numeric (filter out any that failed conversion)
    initial_rows = nrow(df_work)
    df_work = filter(row -> typeof(row.capacity_mwe) <: Number &&
                           typeof(row.occ_usd_per_kw) <: Number, df_work)
    println("Rows after ensuring numeric types: $(nrow(df_work)) (removed $(initial_rows - nrow(df_work)))")

    if nrow(df_work) == 0
        println("\n❌ ERROR: No valid data rows found after cleaning!")
        println("Check that the Excel file has data starting from row 3")
        return nothing
    end

    # Standardize reactor types
    println("\n--- Reactor Type Standardization ---")
    df_work[!, :type] = standardize_reactor_type.(df_work[!, :reactor_type_raw])
    println("Standardized types:")
    println(combine(groupby(df_work, :type), nrow => :count))
    println()

    # Create output dataframe with exact column structure from project_data.csv
    df_output = DataFrame()

    # Column 1: name
    df_output[!, :name] = df_work[!, :name]

    # Column 2: type
    df_output[!, :type] = df_work[!, :type]

    # Column 3: scale (micro/SMR/large)
    if :scale in propertynames(df_work)
        df_output[!, :scale] = df_work[!, :scale]
    else
        # If scale column not found, try to infer from capacity
        println("Warning: scale column not found, inferring from capacity")
        df_output[!, :scale] = map(df_work[!, :capacity_mwe]) do cap
            if ismissing(cap)
                return "Unknown"
            elseif cap < 50
                return "Micro"
            elseif cap <= 300
                return "SMR"
            else
                return "Large"
            end
        end
    end

    # Column 4: country
    if :country in propertynames(df_work)
        df_output[!, :country] = df_work[!, :country]
    else
        println("Warning: country column not found, using 'Unknown'")
        df_output[!, :country] = fill("Unknown", nrow(df_work))
    end

    # Column 5: region (derived from country)
    println("\n--- Regional Categorization ---")
    df_output[!, :region] = country_to_region.(df_output[!, :country])
    println("Regions identified:")
    println(combine(groupby(df_output, :region), nrow => :count))
    println()

    # Column 6: investment (investment cost in USD/MW)
    # Convert OCC from USD/kW to USD/MW (multiply by 1000)
    # Note: Simulation code multiplies this by plant_capacity to get total investment
    df_output[!, :investment] = df_work[!, :occ_usd_per_kw] .* 1000

    # Column 7: plant_capacity (in MWe)
    df_output[!, :plant_capacity] = df_work[!, :capacity_mwe]

    # Column 8: learning_factor (default 0.1 for all)
    df_output[!, :learning_factor] = fill(0.1, nrow(df_work))

    # Column 9: construction_time
    # Use ONLY planned construction time (ignore actual build time)
    println("\n--- Construction Time Debug ---")
    if :construction_time_planned in propertynames(df_work)
        println("construction_time_planned column exists")
        println("Sample values: ", df_work[1:min(5, nrow(df_work)), :construction_time_planned])
        println("Types: ", typeof.(df_work[1:min(5, nrow(df_work)), :construction_time_planned]))
        df_output[!, :construction_time] = coalesce.(df_work[!, :construction_time_planned], 3.0)
    else
        println("construction_time_planned NOT found, using default 3.0")
        df_output[!, :construction_time] = fill(3.0, nrow(df_work))
    end
    println("Final construction_time values: ", df_output[1:min(5, nrow(df_output)), :construction_time])

    # Column 10: operating_time (lifetime in years)
    println("\n--- Operating Time Debug ---")
    if :lifetime_years in propertynames(df_work)
        println("lifetime_years column exists")
        println("Sample values: ", df_work[1:min(5, nrow(df_work)), :lifetime_years])
        println("Types: ", typeof.(df_work[1:min(5, nrow(df_work)), :lifetime_years]))
        df_output[!, :operating_time] = coalesce.(df_work[!, :lifetime_years], 60.0)
    else
        println("lifetime_years column NOT found")
        df_output[!, :operating_time] = fill(60.0, nrow(df_work))
    end
    println("Final operating_time values: ", df_output[1:min(5, nrow(df_output)), :operating_time])

    # Columns 11-12: loadfactor_lower and loadfactor_upper
    # Convert capacity factor from percentage to decimal
    # If capacity factor is provided, use it; otherwise default to 90-95% range
    if :capacity_factor_pct in propertynames(df_work)
        cf_decimal = df_work[!, :capacity_factor_pct] ./ 100.0
        # Use ±2.5% range around the given value, or default to 90-95%
        df_output[!, :loadfactor_lower] = coalesce.(cf_decimal .- 0.025, 0.90)
        df_output[!, :loadfactor_upper] = coalesce.(cf_decimal .+ 0.025, 0.95)
        # Ensure they're in valid range [0, 1]
        df_output[!, :loadfactor_lower] = max.(df_output[!, :loadfactor_lower], 0.0)
        df_output[!, :loadfactor_upper] = min.(df_output[!, :loadfactor_upper], 1.0)
    else
        df_output[!, :loadfactor_lower] = fill(0.90, nrow(df_work))
        df_output[!, :loadfactor_upper] = fill(0.95, nrow(df_work))
    end

    # Column 13: operating_cost_fix (fixed OPEX in USD/MW-yr)
    # Keep as USD/MW-yr - simulation code multiplies by capacity
    if :opex_fixed_usd_per_mw_yr in propertynames(df_work)
        df_output[!, :operating_cost_fix] = coalesce.(df_work[!, :opex_fixed_usd_per_mw_yr], 500.0)
    else
        df_output[!, :operating_cost_fix] = fill(500.0, nrow(df_work))  # Default: 500 USD/MW-yr
    end

    # Column 14: operating_cost_variable (USD/MWh)
    if :opex_variable_usd_per_mwh in propertynames(df_work)
        df_output[!, :operating_cost_variable] = coalesce.(df_work[!, :opex_variable_usd_per_mwh], 2.3326)
    else
        df_output[!, :operating_cost_variable] = fill(2.3326, nrow(df_work))
    end

    # Column 15: operating_cost_fuel (USD/MWh)
    if :fuel_usd_per_mwh in propertynames(df_work)
        df_output[!, :operating_cost_fuel] = coalesce.(df_work[!, :fuel_usd_per_mwh], 6.0)
    else
        df_output[!, :operating_cost_fuel] = fill(6.0, nrow(df_work))
    end

    # Columns 16-17: reference_pj_investment and reference_pj_capacity
    # Based on reactor type
    reference_investment = Float64[]
    reference_capacity = Float64[]

    for rtype in df_output[!, :type]
        ref = get_reference_values(rtype)
        push!(reference_investment, ref.investment)
        push!(reference_capacity, ref.capacity)
    end

    df_output[!, :reference_pj_investment] = reference_investment
    df_output[!, :reference_pj_capacity] = reference_capacity

    # Convert appropriate columns to integers (to match project_data.csv format)
    # These columns should be whole numbers without decimals
    df_output[!, :investment] = round.(Int, df_output[!, :investment])
    df_output[!, :plant_capacity] = round.(Int, df_output[!, :plant_capacity])
    df_output[!, :construction_time] = round.(Int, df_output[!, :construction_time])
    df_output[!, :operating_time] = round.(Int, df_output[!, :operating_time])
    df_output[!, :operating_cost_fix] = round.(Int, df_output[!, :operating_cost_fix])
    df_output[!, :reference_pj_investment] = round.(Int, df_output[!, :reference_pj_investment])
    df_output[!, :reference_pj_capacity] = round.(Int, df_output[!, :reference_pj_capacity])

    # Save as semicolon-delimited CSV
    CSV.write(output_csv, df_output; delim=';')

    println("\n✓ Successfully extracted $(nrow(df_output)) reactors")
    println("✓ Saved to: $output_csv")

    # Print summary statistics
    println("\n--- Summary by Reactor Type ---")
    summary_by_type = combine(groupby(df_output, :type)) do sdf
        (
            count = nrow(sdf),
            avg_capacity = round(mean(sdf.plant_capacity), digits=1),
            avg_investment_per_MW = round(mean(sdf.investment) / 1000, digits=0)  # in thousands USD/MW
        )
    end
    println(summary_by_type)

    # Summary by scale (micro/SMR/large)
    if :scale in propertynames(df_output)
        println("\n--- Summary by Scale (Micro/SMR/Large) ---")
        scale_summary = combine(groupby(df_output, :scale)) do sdf
            (
                count = nrow(sdf),
                avg_capacity = round(mean(sdf.plant_capacity), digits=1),
                avg_investment_per_MW = round(mean(sdf.investment) / 1000, digits=0)  # in thousands USD/MW
            )
        end
        println(scale_summary)
    end

    println("\n--- First 5 reactors ---")
    display_df = df_output[1:min(5, nrow(df_output)), [:name, :type, :scale, :plant_capacity, :investment]]
    println(display_df)

    return df_output
end


# Main execution - runs when script is executed directly from command line
# Usage: julia convert_to_csv.jl
if abspath(PROGRAM_FILE) == @__FILE__
    println("="^60)
    println("Reactor Data Conversion Script")
    println("="^60)
    println()

    excel_file = "_input/reactor_data_raw.xlsx"
    df = extract_reactor_data(excel_file)

    println("\n✓ Conversion complete!")
    println("Output file: _input/reactor_data.csv")
    println("Format matches: project_data.csv structure")
end

# If you're using include() from Julia REPL, call the function directly:
# julia> include("convert_to_csv.jl")
# julia> df = extract_reactor_data("_input/reactor_data_raw.xlsx")
