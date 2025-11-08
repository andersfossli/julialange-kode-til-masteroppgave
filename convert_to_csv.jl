using XLSX
using DataFrames
using CSV

"""
    country_to_region(country)

Map country names to regional categories following Weibezahn et al. (2023)
"""
function country_to_region(country)
    if ismissing(country)
        return "Other"
    end

    country_str = uppercase(strip(string(country)))

    # Emerging Asia: China, India, Pakistan, Russia (checked FIRST)
    emerging_asia = ["CHINA", "RUSSIA", "RUSSIAN FEDERATION", "INDIA", "PAKISTAN"]

    # Western / Developed: USA, Canada, Western Europe, Japan, South Korea
    western_developed = ["USA", "UNITED STATES", "US", "CANADA",
                        "UK", "UNITED KINGDOM", "FRANCE", "FINLAND", "SWEDEN",
                        "GERMANY", "BELGIUM", "NETHERLANDS", "SWITZERLAND",
                        "ITALY", "SPAIN", "JAPAN", "SOUTH KOREA", "KOREA"]

    # Eastern Europe: Former Soviet bloc (excluding Russia)
    eastern_europe = ["UKRAINE", "CZECH REPUBLIC", "SLOVAKIA", "HUNGARY",
                     "POLAND", "BULGARIA", "ROMANIA"]

    # South America
    south_america = ["ARGENTINA", "BRAZIL", "MEXICO", "CHILE"]

    # Middle East / Africa
    middle_east_africa = ["UAE", "SAUDI ARABIA", "IRAN", "TURKEY",
                         "EGYPT", "SOUTH AFRICA"]

    # Check regions in priority order
    if any(occursin(c, country_str) for c in emerging_asia)
        return "Emerging Asia"
    elseif any(occursin(c, country_str) for c in western_developed)
        return "Western / Developed"
    elseif any(occursin(c, country_str) for c in eastern_europe)
        return "Eastern Europe"
    elseif any(occursin(c, country_str) for c in south_america)
        return "South America"
    elseif any(occursin(c, country_str) for c in middle_east_africa)
        return "Middle East / Africa"
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

    pwr_types = ["PWR", "EPR", "AP1000", "APR1400", "VVER1200", "VVER", "IPWR"]
    bwr_types = ["BWR"]
    htr_types = ["HTR", "HTR/GFR", "HTGR"]
    sfr_types = ["SFR", "BN-800", "BN800"]
    msr_types = ["MSR", "MSFR"]
    lfr_types = ["LFR"]
    mr_types = ["MR"]

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
        return string(reactor_type)
    end
end

"""
    get_reference_values(reactor_type)

Get reference project values based on reactor type
"""
function get_reference_values(reactor_type)
    if reactor_type == "PWR"
        return (investment=8600000, capacity=1117)
    elseif reactor_type == "BWR"
        return (investment=9722604, capacity=935)
    elseif reactor_type == "HTR"
        return (investment=7195125, capacity=330)
    elseif reactor_type == "SFR"
        return (investment=5200000, capacity=1250)
    elseif reactor_type == "MSR"
        return (investment=8600000, capacity=1000)
    elseif reactor_type == "LFR"
        return (investment=5200000, capacity=1250)
    else
        return (investment=8600000, capacity=1000)
    end
end

"""
    get_capacity_factor_range(reactor_type)

Get reactor-type-specific capacity factor ranges based on WNA Global Nuclear Industry Performance data.
Uses ±15% ranges from global averages for Monte Carlo simulation uncertainty.

# Capacity Factor Mapping (uniform distribution sampling)
- BWR:     [0.75, 0.95]  # 90% ±15%, capped at 95% (best performing LWR type)
- PWR:     [0.65, 0.95]  # 80% ±15% (most common LWR type)
- HTR:     [0.55, 0.85]  # 70% ±15% (gas-cooled reactor, limited operational experience)
- HTGR:    [0.65, 0.95]  # 80% ±15% (using default - overlaps with HTR category)
- HTR/GFR: [0.55, 0.85]  # 70% ±15% (gas-cooled fast reactor, no operational data)
- SFR:     [0.55, 0.85]  # 70% ±15% (fast breeder, limited operational experience)
- LFR:     [0.55, 0.85]  # 70% ±15% (fast breeder category, no operational data)
- MSR:     [0.65, 0.95]  # 80% ±15% (default - no operational data)
- MSFR:    [0.65, 0.95]  # 80% ±15% (default - no operational data)
- MR:      [0.65, 0.95]  # 80% ±15% (default - no operational data)

# Returns
- Tuple (CF_min, CF_max) for uniform distribution sampling
"""
function get_capacity_factor_range(reactor_type)
    # Reactor-type-specific capacity factor ranges from WNA data
    cf_ranges = Dict(
        "BWR"     => (0.75, 0.95),  # 90% ±15%, capped at 95%
        "PWR"     => (0.65, 0.95),  # 80% ±15%
        "HTR"     => (0.55, 0.85),  # 70% ±15% (gas-cooled)
        "HTGR"    => (0.65, 0.95),  # 80% ±15% (using default)
        "HTR/GFR" => (0.55, 0.85),  # 70% ±15% (gas-cooled fast)
        "SFR"     => (0.55, 0.85),  # 70% ±15% (fast breeder)
        "LFR"     => (0.55, 0.85),  # 70% ±15% (fast breeder)
        "MSR"     => (0.65, 0.95),  # 80% ±15% (default)
        "MSFR"    => (0.65, 0.95),  # 80% ±15% (default)
        "MR"      => (0.65, 0.95)   # 80% ±15% (default)
    )

    # Error handling: ensure reactor type is mapped
    if !haskey(cf_ranges, reactor_type)
        error("Unmapped reactor type: '$reactor_type'. " *
              "Valid types: $(join(keys(cf_ranges), ", "))")
    end

    return cf_ranges[reactor_type]
end

"""
    clean_numeric_string(s)

Clean numeric strings by removing spaces, dashes, and 'x' markers
Handle ranges like "85-90" by taking midpoint
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

Extract reactor data from Excel and convert to CSV format
"""
function extract_reactor_data(excel_file::String; output_csv::String="_input/reactor_data.csv")
    println("Reading Excel file: $excel_file")
    
    # Read Excel
    xf = XLSX.readxlsx(excel_file)
    sheet = xf["Datasheet"]
    df = DataFrame(XLSX.gettable(sheet; first_row=3))

    # Column mapping
    column_mapping = Dict(
        "Project" => :name,
        "Type" => :reactor_type_raw,
        "Country" => :country,
        "Capacity (net MWe)" => :capacity_mwe,
        "OCC (USD2020/kW)" => :occ_usd_per_kw,
        "planned Construction time (y)" => :construction_time_planned,
        "Actual / total build time (y)" => :construction_time_actual,
        "Lifetime (y)" => :lifetime_years,
        "Capacity factor (%)" => :capacity_factor_pct,
        "OPEX fixed (USD2020/MW-yr)" => :opex_fixed_usd_per_mw_yr,
        "OPEX variable (USD2020/MWh)" => :opex_variable_usd_per_mwh,
        "Fuel (USD2020/MWh)" => :fuel_usd_per_mwh,
        "Scale" => :scale,
        "scale" => :scale,
        "Large medium micro" => :scale,
        "FOAK/NOAK*" => :foak_noak,
        "Year" => :year
    )

    # Select and rename available columns
    available_cols = Dict{String, Symbol}()
    for (old_name, new_name) in column_mapping
        if old_name in names(df)
            available_cols[old_name] = new_name
        end
    end

    df_work = select(df, keys(available_cols)...)
    rename!(df_work, available_cols)
    rename!(df_work, Symbol.(names(df_work)))

    # Filter out header rows
    if :name in propertynames(df_work)
        df_work = filter(row -> !ismissing(row.name) &&
                                !occursin(r"^Project$"i, string(row.name)), df_work)
    end
    if :reactor_type_raw in propertynames(df_work)
        df_work = filter(row -> !ismissing(row.reactor_type_raw) &&
                                !occursin(r"^Type$"i, string(row.reactor_type_raw)), df_work)
    end

    # Clean numeric columns
    numeric_cols = [:capacity_mwe, :occ_usd_per_kw, :construction_time_planned,
                    :construction_time_actual, :lifetime_years, :capacity_factor_pct,
                    :opex_fixed_usd_per_mw_yr, :opex_variable_usd_per_mwh, :fuel_usd_per_mwh,
                    :year]

    for col in numeric_cols
        if col in propertynames(df_work)
            df_work[!, col] = clean_numeric_string.(df_work[!, col])
        end
    end

    # Remove rows with missing critical data
    df_work = dropmissing(df_work, [:name, :capacity_mwe, :occ_usd_per_kw])

    # Ensure numeric columns are numeric
    df_work = filter(row -> typeof(row.capacity_mwe) <: Number &&
                           typeof(row.occ_usd_per_kw) <: Number, df_work)

    if nrow(df_work) == 0
        error("No valid data rows found after cleaning!")
    end

    # Standardize reactor types
    df_work[!, :type] = standardize_reactor_type.(df_work[!, :reactor_type_raw])

    # Remove reactors without operational reference
    reactors_to_exclude = ["IMSR (300)", "SSR-W", "e-Vinci", "BREST-OD-300", "Brest-OD-300"]
    df_work = filter(row -> !(row.name in reactors_to_exclude), df_work)

    # Create output dataframe
    df_output = DataFrame()

    # Column 1: name
    df_output[!, :name] = df_work[!, :name]

    # Column 2: type
    df_output[!, :type] = df_work[!, :type]

    # Column 3: scale
    if :scale in propertynames(df_work)
        df_output[!, :scale] = df_work[!, :scale]
    else
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
        df_output[!, :country] = fill("Unknown", nrow(df_work))
    end

    # Column 5: region
    df_output[!, :region] = country_to_region.(df_output[!, :country])

    # Column 6: investment (USD/MW)
    df_output[!, :investment] = df_work[!, :occ_usd_per_kw] .* 1000

    # Column 7: plant_capacity (MWe)
    df_output[!, :plant_capacity] = df_work[!, :capacity_mwe]

    # Column 8: learning_factor
    df_output[!, :learning_factor] = fill(0.1, nrow(df_work))

    # Column 9: construction_time
    if :construction_time_planned in propertynames(df_work)
        df_output[!, :construction_time] = coalesce.(df_work[!, :construction_time_planned], 3.0)
    else
        df_output[!, :construction_time] = fill(3.0, nrow(df_work))
    end

    # Column 10: operating_time
    if :lifetime_years in propertynames(df_work)
        df_output[!, :operating_time] = coalesce.(df_work[!, :lifetime_years], 60.0)
    else
        df_output[!, :operating_time] = fill(60.0, nrow(df_work))
    end

    # Columns 11-12: loadfactor range (reactor-type-specific from WNA data)
    # Use reactor-type-specific capacity factor ranges based on operational performance
    loadfactor_lower = Float64[]
    loadfactor_upper = Float64[]

    for rtype in df_output[!, :type]
        cf_min, cf_max = get_capacity_factor_range(rtype)
        push!(loadfactor_lower, cf_min)
        push!(loadfactor_upper, cf_max)
    end

    df_output[!, :loadfactor_lower] = loadfactor_lower
    df_output[!, :loadfactor_upper] = loadfactor_upper

    # Column 13: operating_cost_fix
    if :opex_fixed_usd_per_mw_yr in propertynames(df_work)
        df_output[!, :operating_cost_fix] = coalesce.(df_work[!, :opex_fixed_usd_per_mw_yr], 500.0)
    else
        df_output[!, :operating_cost_fix] = fill(500.0, nrow(df_work))
    end

    # Column 14: operating_cost_variable
    if :opex_variable_usd_per_mwh in propertynames(df_work)
        df_output[!, :operating_cost_variable] = coalesce.(df_work[!, :opex_variable_usd_per_mwh], 2.3326)
    else
        df_output[!, :operating_cost_variable] = fill(2.3326, nrow(df_work))
    end

    # Column 15: operating_cost_fuel
    if :fuel_usd_per_mwh in propertynames(df_work)
        df_output[!, :operating_cost_fuel] = coalesce.(df_work[!, :fuel_usd_per_mwh], 6.0)
    else
        df_output[!, :operating_cost_fuel] = fill(6.0, nrow(df_work))
    end

    # Columns 16-17: reference values
    reference_investment = Float64[]
    reference_capacity = Float64[]

    for rtype in df_output[!, :type]
        ref = get_reference_values(rtype)
        push!(reference_investment, ref.investment)
        push!(reference_capacity, ref.capacity)
    end

    df_output[!, :reference_pj_investment] = reference_investment
    df_output[!, :reference_pj_capacity] = reference_capacity

    # Column 18: year
    if :year in propertynames(df_work)
        df_output[!, :year] = coalesce.(df_work[!, :year], 2020)
    else
        df_output[!, :year] = fill(2020, nrow(df_work))
    end

    # Convert to integers where appropriate
    df_output[!, :investment] = round.(Int, df_output[!, :investment])
    df_output[!, :plant_capacity] = round.(Int, df_output[!, :plant_capacity])
    df_output[!, :construction_time] = round.(Int, df_output[!, :construction_time])
    df_output[!, :operating_time] = round.(Int, df_output[!, :operating_time])
    df_output[!, :operating_cost_fix] = round.(Int, df_output[!, :operating_cost_fix])
    df_output[!, :reference_pj_investment] = round.(Int, df_output[!, :reference_pj_investment])
    df_output[!, :reference_pj_capacity] = round.(Int, df_output[!, :reference_pj_capacity])
    df_output[!, :year] = round.(Int, df_output[!, :year])

    # Save CSV
    CSV.write(output_csv, df_output; delim=';')

    println("\n✓ Successfully converted $(nrow(df_output)) reactors")
    println("✓ Saved to: $output_csv")

    return df_output
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("="^60)
    println("Reactor Data Conversion: Excel → CSV")
    println("="^60)

    excel_file = "_input/reactor_data_raw.xlsx"
    df = extract_reactor_data(excel_file)

    println("\n✓ Conversion complete!")
end