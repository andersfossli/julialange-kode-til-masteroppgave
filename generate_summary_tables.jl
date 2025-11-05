using CSV
using DataFrames
using Statistics

# Try to load PrettyTables for nice formatting
const PRETTYTABLES_AVAILABLE = try
    using PrettyTables
    true
catch
    @warn "PrettyTables not installed. Install with: import Pkg; Pkg.add(\"PrettyTables\")"
    false
end

"""
    generate_summary_tables(csv_file::String="_input/reactor_data.csv")

Read reactor data CSV and generate summary statistics tables
Saves tables to _output/ directory
"""
function generate_summary_tables(csv_file::String="_input/reactor_data.csv")
    
    if !isfile(csv_file)
        error("CSV file not found: $csv_file")
    end

    # Read CSV
    df = CSV.File(csv_file, delim=';') |> DataFrame
    total_reactors = nrow(df)

    println("\n" * "="^80)
    println("REACTOR DATASET SUMMARY STATISTICS")
    println("="^80)
    println("Total reactors: $total_reactors")
    println("Source: $csv_file")
    println()

    # ========== TABLE 1: BY SCALE ==========
    scale_data = DataFrame(
        Category = String[],
        Range = String[],
        Count = Int[],
        Share = String[],
        Description = String[]
    )

    scale_order = ["Large", "SMR", "Micro"]
    scale_ranges = Dict(
        "Large" => ">300 MWe",
        "SMR" => "50-300 MWe",
        "Micro" => "<50 MWe"
    )
    scale_descriptions = Dict(
        "Large" => "Conventional Gen III/III+ units",
        "SMR" => "Modular designs under development",
        "Micro" => "Small-scale/off-grid designs"
    )

    scale_counts = combine(groupby(df, :scale), nrow => :count)
    for scale in scale_order
        row = filter(r -> r.scale == scale, scale_counts)
        if nrow(row) > 0
            count = row[1, :count]
            share = "$(round(count / total_reactors * 100, digits=1))%"
            push!(scale_data, (scale, scale_ranges[scale], count, share, scale_descriptions[scale]))
        end
    end

    println("="^80)
    println("TABLE 1: REACTORS BY SCALE")
    println("="^80)
    if PRETTYTABLES_AVAILABLE
        pretty_table(scale_data,
                     backend=Val(:text),
                     header=["Category", "Range", "Count", "Share", "Description"],
                     alignment=[:l, :l, :r, :r, :l],
                     crop=:none)
    else
        println(scale_data)
    end

    # ========== TABLE 2: BY TYPE ==========
    type_data = DataFrame(
        Type = String[],
        Full_Name = String[],
        Count = Int[],
        Share = String[],
        Description = String[]
    )

    type_full_names = Dict(
        "PWR" => "Pressurized Water Reactor",
        "BWR" => "Boiling Water Reactor",
        "SFR" => "Sodium-cooled Fast Reactor",
        "HTR" => "High-Temperature Gas-cooled Reactor"
    )

    type_descriptions = Dict(
        "PWR" => "Dominant type in SMR and Large",
        "BWR" => "GE, Hitachi designs",
        "SFR" => "BN, CEFR, ARC designs",
        "HTR" => "Fort St. Vrain, HTR-PM, EM2"
    )

    type_counts = combine(groupby(df, :type), nrow => :count)
    sort!(type_counts, :count, rev=true)

    for row in eachrow(type_counts)
        rtype = row.type
        count = row.count
        share = "$(round(count / total_reactors * 100, digits=1))%"
        full_name = get(type_full_names, rtype, rtype)
        desc = get(type_descriptions, rtype, "")
        push!(type_data, (rtype, full_name, count, share, desc))
    end

    println("\n" * "="^80)
    println("TABLE 2: REACTORS BY TYPE")
    println("="^80)
    if PRETTYTABLES_AVAILABLE
        pretty_table(type_data,
                     backend=Val(:text),
                     header=["Type", "Full Name", "Count", "Share", "Description"],
                     alignment=[:l, :l, :r, :r, :l],
                     crop=:none)
    else
        println(type_data)
    end

    # ========== TABLE 3: BY REGION ==========
    region_data = DataFrame(
        Region = String[],
        Count = Int[],
        Share = String[],
        Countries = String[],
        Description = String[]
    )

    region_countries = Dict(
        "Western / Developed" => "USA, UK, France, Canada, Japan, S. Korea",
        "Emerging Asia" => "China, Russia, India, Pakistan",
        "South America" => "Argentina, Brazil",
        "Middle East / Africa" => "UAE, Saudi Arabia, Egypt, S. Africa"
    )

    region_descriptions = Dict(
        "Western / Developed" => "OECD vendor-led projects",
        "Emerging Asia" => "State-driven programs",
        "South America" => "Regional development",
        "Middle East / Africa" => "Emerging markets"
    )

    region_counts = combine(groupby(df, :region), nrow => :count)
    sort!(region_counts, :count, rev=true)

    for row in eachrow(region_counts)
        region = row.region
        count = row.count
        share = "$(round(count / total_reactors * 100, digits=1))%"
        countries = get(region_countries, region, "Various")
        desc = get(region_descriptions, region, "")
        push!(region_data, (region, count, share, countries, desc))
    end

    println("\n" * "="^80)
    println("TABLE 3: REACTORS BY GEOGRAPHIC REGION")
    println("="^80)
    if PRETTYTABLES_AVAILABLE
        pretty_table(region_data,
                     backend=Val(:text),
                     header=["Region", "Count", "Share", "Countries", "Description"],
                     alignment=[:l, :r, :r, :l, :l],
                     crop=:none)
    else
        println(region_data)
    end

    println("\n" * "="^80)

    # Save tables to CSV
    CSV.write("_output/summary_by_scale.csv", scale_data)
    CSV.write("_output/summary_by_type.csv", type_data)
    CSV.write("_output/summary_by_region.csv", region_data)

    println("\n✓ Summary tables saved to _output/")
    println("  - summary_by_scale.csv")
    println("  - summary_by_type.csv")
    println("  - summary_by_region.csv")

    return (scale=scale_data, type=type_data, region=region_data)
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("="^60)
    println("Reactor Dataset Summary Tables")
    println("="^60)

    tables = generate_summary_tables()

    println("\n✓ Table generation complete!")
end