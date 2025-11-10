"""
The investment_plot function generates a comparison plot of investment values for manufacturer vs. estimated costs for multiple projects. The function takes in two arguments: a vector pjs containing project objects to be compared, and a boolean scaling_plot which indicates whether to scale the investment values by the project capacity.
The function first generates scaled investment values for all projects using the gen_scaled_investment function, and stores them in a DataFrame. The DataFrame columns are named after the project names.
The function then creates a new figure object and defines the x-axis label and y-axis tick labels based on the input DataFrame.
The function then plots the investment values using three different plot types: scatter for the manufacturer's investment values, and rangebars for the other two project's investment values. The rangebars represent the minimum and maximum investment values for each project.
Finally, the function adds a legend to the plot and returns the figure object.
"""
# DEAD CODE: Commented out - replaced by investment_plot_by_scale()
# This function showed all 15 reactors in one plot (too crowded)
# The new function groups by scale (Micro/SMR/Large) for better visualization
# function investment_plot(pjs, scaling_plot)
#     scaled_investments = DataFrame();
#     for p in eachindex(pjs)
#         scaled_investments.res = gen_scaled_investment(scaling_plot, pjs[p])
#         rename!(scaled_investments,:res => pjs[p].name)
#     end
#     xlabel = "[USD/MW]";
#     yticks = names(scaled_investments);
#     fig_invest_comparison = Figure();
#     ax_invest = Axis(fig_invest_comparison[1,1], yticks = (1:length(yticks), yticks), xscale = log10, xlabel = xlabel);
#     rothwell = rangebars!(ax_invest, 1.2:length(yticks)+0.2, collect(scaled_investments[5,:]), collect(scaled_investments[4,:]), linewidth = 6, whiskerwidth = 12, direction = :x, transparency = :true, color = :green)
#     roulstone = rangebars!(ax_invest, 1:length(yticks), collect(scaled_investments[3,:]), collect(scaled_investments[2,:]), linewidth = 6, whiskerwidth = 12, direction = :x, transparency = :true, color = :darkblue)
#     manufacturer = scatter!(ax_invest, collect(scaled_investments[1,:]), 1:length(yticks), marker = :star5, color = :red)
#     Legend(fig_invest_comparison[1, 1],
#         [manufacturer, roulstone, rothwell],
#         ["Manufacturer", "Roulstone", "Rothwell"],
#         tellheight = false,
#         tellwidth = false,
#         halign = :right, valign = :bottom,
#         framevisible = false, orientation = :vertical)
#     return fig_invest_comparison
# end

"""
The hist_invest_plot function generates a histogram plot of the investment values of a given project for a specified number of trials. The function takes in six arguments:
- `n::Int64`: The number of trials.
- `wacc::Vector`: A vector of weighted average cost of capital values for each trial (not used for the plot but needed for the gen_rand_vars function).
- `electricity_price_mean::Float64`: Mean electricity price (not used for the plot but needed for the gen_rand_vars function).
- `pj::project`: An instance of a project class that contains information about the project.
- `i::Int64`: The row index of the subplot where the histogram plot will be added.
- `j::Int64`: The column index of the subplot where the histogram plot will be added.
- `hist_invest = Figure()`: The histogram plot is added to this figure object. If not provided, a new figure object will be created.

The function first generates random variables using the gen_rand_vars function for two different options of scaling. It then adds the histograms to the plot. Each histogram is normalized to show the probability density function.
"""
function hist_invest_plot(n::Int64, wacc::Vector, electricity_price_mean::Float64, pj::project, i::Int64, j::Int64, hist_invest = Figure())

    random_vars_roulstone = gen_rand_vars(opts_scaling[2], n, wacc, electricity_price_mean, pj)
    random_vars_rothwell = gen_rand_vars(opts_scaling[3], n, wacc, electricity_price_mean, pj)
    # random_vars_uniform = gen_rand_vars(opts_scaling[4], n, wacc, electricity_price_mean, pj)

    ax_invest = Axis(hist_invest[i,j]);
    hidedecorations!(ax_invest)
    hist!(ax_invest, vec(random_vars_rothwell.investment), normalization = :probability, color = (:green, 0.8), strokewidth = 1, strokecolor = :black)
    hist!(ax_invest, vec(random_vars_roulstone.investment), normalization = :probability, color = (:darkblue, 0.8), strokewidth = 1, strokecolor = :black)
    # hist!(ax_invest, vec(random_vars_uniform.investment), normalization = :probability, color = (:gray, 0.5), strokewidth = 1, strokecolor = :black)
    Label(hist_invest[i, j, Top()], pj.name, font = "Noto Sans Bold", padding = (0, 6, 6, 0))
    Label(hist_invest[i, j], pj.type, font = "Noto Sans Bold", fontsize = 8, halign = :right, valign = :top, justification = :center, padding = (4, 6, 6, 4), tellheight = false, tellwidth = false)
    
    return hist_invest

end

"""
The density_invest_plot function generates a probability density plot of the investment values of a given project for a specified number of trials. The function takes in six arguments:
- `n::Int64`: The number of trials.
- `wacc::Vector`: A vector of weighted average cost of capital values for each trial (not used for the plot but needed for the gen_rand_vars function).
- `electricity_price_mean::Float64`: Mean electricity price (not used for the plot but needed for the gen_rand_vars function).
- `pj::project`: An instance of a project class that contains information about the project.
- `i::Int64`: The row index of the subplot where the histogram plot will be added.
- `j::Int64`: The column index of the subplot where the histogram plot will be added.
- `density_invest = Figure()`: The histogram plot is added to this figure object. If not provided, a new figure object will be created.

The function first generates random variables using the gen_rand_vars function for two different options of scaling. It then adds the probability density function to the plot.
"""
function density_invest_plot(n::Int64, wacc::Vector, electricity_price_mean::Float64, pj::project, i::Int64, j::Int64, density_invest = Figure())

    random_vars_roulstone = gen_rand_vars(opts_scaling[2], n, wacc, electricity_price_mean, pj)
    random_vars_rothwell = gen_rand_vars(opts_scaling[3], n, wacc, electricity_price_mean, pj)
    # random_vars_uniform = gen_rand_vars(opts_scaling[4], n, wacc, electricity_price_mean, pj)

    ax_invest = Axis(density_invest[i,j]);
    hidedecorations!(ax_invest)
    density!(ax_invest, vec(random_vars_rothwell.investment), color = (:green, 0.8), strokewidth = 1, strokecolor = :black)
    density!(ax_invest, vec(random_vars_roulstone.investment), color = (:darkblue, 0.8), strokewidth = 1, strokecolor = :black)
    # density!(ax_invest, vec(random_vars_uniform.investment), color = (:gray, 0.5), strokewidth = 1, strokecolor = :black)
    Label(density_invest[i, j, Top()], pj.name, font = "Noto Sans Bold", padding = (0, 6, 6, 0))
    Label(density_invest[i, j], pj.type, font = "Noto Sans Bold", fontsize = 8, halign = :right, valign = :top, justification = :center, padding = (4, 6, 6, 4), tellheight = false, tellwidth = false)

    return density_invest
    
end

"""
The mcs_plot function creates a box plot figure for Monte Carlo simulation results organized by three categories of nuclear reactor types: BWR & PWR, HTR, and SFR. The input arguments are:
    mcs_results: a DataFrame containing the Monte Carlo simulation results. The rows represent individual runs, and the columns represent different reactor types and parameters.
    title: a string specifying the title of the figure.
    ylabel: a string specifying the label of the y-axis.
    pjs: (optional) Vector of reactor project objects to dynamically group by type and scale
The function generates a box plot for each of the three categories of reactor types. The x-axis of each box plot is labeled with the parameter names, and the y-axis is labeled with the ylabel argument. The color of the box plots depends on the reactor type. The figure has a title specified by the title argument, and each box plot has a label specified by the corresponding reactor type.
"""
function mcs_plot(mcs_results, title::String, ylabel::String, pjs::Union{Vector,Nothing}=nothing)

    # If pjs is provided, dynamically group reactors by type
    if !isnothing(pjs)
        # Group reactor names by type
        pwr_bwr_reactors = String[]
        htr_reactors = String[]
        sfr_reactors = String[]

        for pj in pjs
            if pj.name in names(mcs_results)
                if pj.type == "PWR" || pj.type == "BWR"
                    push!(pwr_bwr_reactors, pj.name)
                elseif pj.type == "HTR"
                    push!(htr_reactors, pj.name)
                elseif pj.type == "SFR"
                    push!(sfr_reactors, pj.name)
                end
            end
        end

        xticks_wc = pwr_bwr_reactors
        xticks_ht = htr_reactors
        xticks_sf = sfr_reactors
    else
        # Fallback to old hardcoded behavior
        xticks_wc = names(mcs_results)[1:9];
        xticks_ht = names(mcs_results)[10:12];
        xticks_sf = names(mcs_results)[13:15];
    end

    mcs_boxplot = Figure();
    n = nrow(mcs_results)  # Number of Monte Carlo samples

    # Calculate 10-90% quantile range for y-axis limits (across all reactors)
    all_data = Float64[]
    for reactor_name in names(mcs_results)
        if reactor_name in vcat(xticks_wc, xticks_ht, xticks_sf)
            append!(all_data, mcs_results[!, reactor_name])
        end
    end
    q10 = quantile(all_data, 0.10)
    q90 = quantile(all_data, 0.90)
    range_width = q90 - q10
    ylim_low = q10 - 0.1 * range_width
    ylim_high = q90 + 0.1 * range_width

    ax_wc = Axis(mcs_boxplot[1,1],
                 xticks = (1:length(xticks_wc), xticks_wc),
                 ylabel = ylabel,
                 limits = (nothing, nothing, ylim_low, ylim_high));
    ax_wc.xticklabelrotation = π / 3;
    ax_wc.yticklabelrotation = π / 2;
    ax_wc.xticklabelalign = (:right, :center);

    for (i, reactor_name) in enumerate(xticks_wc)
        boxplot!(ax_wc, fill(i,n), mcs_results[!,reactor_name], color = :orangered)
    end;

    Label(mcs_boxplot[1, 1, Top()], "BWR & PWR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0))

    ax_ht = Axis(mcs_boxplot[1,2],
                 xticks = (1:length(xticks_ht), xticks_ht),
                 limits = (nothing, nothing, ylim_low, ylim_high));
    ax_ht.xticklabelrotation = π / 3;
    ax_ht.yticklabelrotation = π / 2;
    ax_ht.xticklabelalign = (:right, :center);

    for (i, reactor_name) in enumerate(xticks_ht)
        boxplot!(ax_ht, fill(i,n), mcs_results[!,reactor_name], color = :gold)
    end;

    Label(mcs_boxplot[1, 2, Top()], "HTR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    ax_sf = Axis(mcs_boxplot[1,3],
                 xticks = (1:length(xticks_sf), xticks_sf),
                 limits = (nothing, nothing, ylim_low, ylim_high));
    ax_sf.xticklabelrotation = π / 3;
    ax_sf.yticklabelrotation = π / 2;
    ax_sf.xticklabelalign = (:right, :center);

    for (i, reactor_name) in enumerate(xticks_sf)
        boxplot!(ax_sf, fill(i,n), mcs_results[!,reactor_name], color = :teal)
    end;

    Label(mcs_boxplot[1, 3, Top()], "SFR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    Label(mcs_boxplot[0, :], title, fontsize = 24, font = "Noto Sans Bold", color = (:black, 0.25))

    # Adjust column sizes based on number of reactors in each category
    total_reactors = length(xticks_wc) + length(xticks_ht) + length(xticks_sf)
    colsize!(mcs_boxplot.layout, 1, Relative(length(xticks_wc) / total_reactors));
    colsize!(mcs_boxplot.layout, 2, Relative(length(xticks_ht) / total_reactors));
    colsize!(mcs_boxplot.layout, 3, Relative(length(xticks_sf) / total_reactors));

    return mcs_boxplot

end

"""
The function si_plot() takes in two arguments, si_results and title, and returns a heatmap figure displaying sensitivity indices.
The si_results argument is expected to be a DataFrame containing the results of a sensitivity analysis, where each row corresponds to a different input variable and each column contains sensitivity indices for that variable.
The title argument is a string used as the title of the resulting heatmap figure.
The function filters the si_results DataFrame into two separate DataFrames, one for first-order sensitivity indices (si_s) and one for total-order sensitivity indices (si_st), and extracts the variable names and index values for the first variable column.
The function then creates a new figure si_heatmap and initializes three axes. The first axis (ax_s) displays the first-order sensitivity indices heatmap, the second axis displays the colorbar for the first-order sensitivity indices, and the third axis (ax_st) displays the total-order sensitivity indices heatmap.
The function populates the first-order and total-order sensitivity indices heatmaps with data using the heatmap! function from the Plots package, sets the tick labels and rotation, and adds text annotations showing the rounded sensitivity index values.
The function also adds labels to the first-order and total-order sensitivity indices heatmaps and the title to the top of the figure. Finally, the function returns the resulting si_heatmap figure.
"""

function si_plot(si_results, title::String)

    si_s = filter(:si => index -> index == "S", si_results);
    si_st = filter(:si => index -> index == "ST", si_results);
    xticks = si_s.var;
    yticks = names(si_results)[3:end];
    data_s = Matrix(si_s[:,3:end]);
    data_st = Matrix(si_st[:,3:end]);

    si_heatmap = Figure();

    ax_s = Axis(si_heatmap[1, 1], xticks = (1:length(xticks), xticks), yticks = (1:length(yticks), yticks));
    ax_s.xticklabelrotation = π / 3;
    ax_s.xticklabelalign = (:right, :center);
    hmap_s = heatmap!(ax_s, data_s, colormap = :deep, colorrange = (0, 1));
    for i in 1:length(xticks), j in 1:length(yticks)
        txtcolor = data_s[i, j] > 0.5 ? :white : :black
        text!(ax_s, "$(round(data_s[i,j], digits = 2))", position = (i, j),
            color = txtcolor, fontsize = 12, align = (:center, :center))
    end
    Label(si_heatmap[1, 1, Top()], "first-order effect", font = "Noto Sans Bold", padding = (0, 6, 6, 0))

    Colorbar(si_heatmap[1, 2], hmap_s; width = 15, ticksize = length(yticks));

    ax_st = Axis(si_heatmap[1, 3], xticks = (1:length(xticks), xticks), yticks = (1:length(yticks), yticks), yaxisposition = :right);
    ax_st.xticklabelrotation = π / 3;
    ax_st.xticklabelalign = (:right, :center);
    hmap_st = heatmap!(ax_st, data_st, colormap = :deep, colorrange = (0, 1));
    for i in 1:length(xticks), j in 1:length(yticks)
        txtcolor = data_st[i, j] > 0.5 ? :white : :black
        text!(ax_st, "$(round(data_st[i,j], digits = 2))", position = (i, j),
            color = txtcolor, fontsize = 12, align = (:center, :center))
    end
    Label(si_heatmap[1, 3, Top()], "total-order effect", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    Label(si_heatmap[0, :], title, fontsize = 24, font = "Noto Sans Bold", color = (:black, 0.25));

    return si_heatmap

end
"""
The lcoe_scale_histogram function creates three separate histograms showing LCOE distributions
by reactor scale (Micro/SMR/Large) with mean and median lines for each group.

Arguments:
- `lcoe_results::DataFrame`: DataFrame containing LCOE results from Monte Carlo simulation (one column per reactor)
- `pjs::Vector`: Vector of project objects containing reactor information including scale

Returns:
- Figure object with three subplots (one per scale)
"""
function lcoe_scale_histogram(lcoe_results::DataFrame, pjs::Vector)

    # Group reactors by scale
    scale_groups = Dict("Micro" => Int[], "SMR" => Int[], "Large" => Int[])

    for (idx, pj) in enumerate(pjs)
        if haskey(scale_groups, pj.scale)
            push!(scale_groups[pj.scale], idx)
        end
    end

    # Collect LCOE data for each scale group
    lcoe_by_scale = Dict{String, Vector{Float64}}()

    for (scale, indices) in scale_groups
        if !isempty(indices)
            # Concatenate all LCOE values for this scale
            scale_data = Float64[]
            for idx in indices
                reactor_name = pjs[idx].name
                if reactor_name in names(lcoe_results)
                    append!(scale_data, lcoe_results[!, reactor_name])
                end
            end
            lcoe_by_scale[scale] = scale_data
        end
    end

    # Create figure with 3 subplots (1 row, 3 columns)
    fig = Figure(size = (1800, 500))

    # Define colors for each scale
    colors = Dict("Micro" => :red, "SMR" => :blue, "Large" => :green)

    # Define scale order and column positions
    scale_order = ["Micro", "SMR", "Large"]

    for (col_idx, scale) in enumerate(scale_order)
        if haskey(lcoe_by_scale, scale) && !isempty(lcoe_by_scale[scale])
            data = lcoe_by_scale[scale]

            # Calculate 10th and 90th percentiles for x-axis limits (removes extreme tails)
            q10 = quantile(data, 0.10)
            q90 = quantile(data, 0.90)

            # Add 10% padding to the range for better visualization
            range_width = q90 - q10
            xlim_low = q10 - 0.1 * range_width
            xlim_high = q90 + 0.1 * range_width

            # Create axis for this subplot with limited x-range
            ax = Axis(fig[1, col_idx],
                      xlabel = "LCOE [USD/MWh]",
                      ylabel = "Probability Density",
                      title = "$scale Reactors (10-90% quantile range)",
                      limits = (xlim_low, xlim_high, nothing, nothing))

            # Plot histogram with probability density normalization
            hist!(ax, data,
                  bins = 50,
                  normalization = :pdf,  # Probability density function
                  color = (colors[scale], 0.6),
                  strokewidth = 1,
                  strokecolor = colors[scale])

            # Calculate and plot mean
            mean_val = mean(data)
            vlines!(ax, [mean_val],
                    color = :black,
                    linestyle = :solid,
                    linewidth = 2,
                    label = "Mean: $(round(mean_val, digits=1))")

            # Calculate and plot median
            median_val = median(data)
            vlines!(ax, [median_val],
                    color = :black,
                    linestyle = :dash,
                    linewidth = 2,
                    label = "Median: $(round(median_val, digits=1))")

            # Add legend
            axislegend(ax, position = :rt)
        end
    end

    # Add overall title
    Label(fig[0, :], "LCOE Distribution by Reactor Scale",
          fontsize = 20, font = "Noto Sans Bold")

    return fig
end

"""
    learning_curve_plot(outputpath::String, opt_scaling::String, learning_scenarios::Vector;
                        reactor_name::Union{Nothing,String}=nothing, scale_filter::Union{Nothing,String}=nothing)

Create a plot showing mean LCOE vs. number of units (N) for different learning scenarios.

# Arguments
- `outputpath::String`: Path to directory containing CSV output files
- `opt_scaling::String`: Scaling method used (e.g., "roulstone", "manufacturer")
- `learning_scenarios::Vector`: Vector of tuples (N_unit, tag) matching the scenarios run
- `reactor_name::Union{Nothing,String}=nothing`: Optional specific reactor to plot (default: average across all)
- `scale_filter::Union{Nothing,String}=nothing`: Optional filter by scale ("Micro", "SMR", "Large")

# Returns
- Figure object with learning curve plot

# Example
```julia
learning_scenarios = [(1, "LR10_N1_k120"), (2, "LR10_N2_k120"), (4, "LR10_N4_k120"), (6, "LR10_N6_k120")]
fig = learning_curve_plot("_output", "roulstone", learning_scenarios)
save("_output/fig-learning_curve.pdf", fig)
```
"""
function learning_curve_plot(outputpath::String, opt_scaling::String, learning_scenarios::Vector;
                             reactor_name::Union{Nothing,String}=nothing,
                             scale_filter::Union{Nothing,String}=nothing)

    # Collect mean LCOE for each scenario
    N_values = Int[]
    mean_lcoe = Float64[]
    scenario_labels = String[]

    for (N_unit, tag) in learning_scenarios
        # Read LCOE summary for this scenario
        summary_file = "$outputpath/mcs-lcoe_summary-$opt_scaling-$tag.csv"

        if !isfile(summary_file)
            @warn("File not found: $summary_file, skipping scenario N=$N_unit")
            continue
        end

        df_summary = CSV.File(summary_file) |> DataFrame

        # Calculate mean LCOE across reactors (or for specific reactor)
        if isnothing(reactor_name)
            # Average across all reactors
            # In describe() output: each row is a reactor, "mean" is a column
            if !("mean" in names(df_summary))
                @warn("'mean' column not found in $summary_file, skipping scenario N=$N_unit")
                continue
            end

            # Apply scale filter if specified
            if !isnothing(scale_filter)
                # Filter rows by reactor names that match the scale
                # This is a simplified approach - you may need reactor metadata
                df_summary = filter(row -> occursin(scale_filter, string(row.variable)), df_summary)
            end

            # Get mean values for all reactors (skip missing values)
            mean_values = filter(!ismissing, df_summary[!, :mean])
            if isempty(mean_values)
                @warn("No valid mean values found in $summary_file")
                continue
            end
            lcoe_mean = mean(mean_values)
        else
            # Specific reactor - find the row for this reactor
            reactor_row = filter(row -> row.variable == reactor_name, df_summary)
            if nrow(reactor_row) == 0
                @warn("Reactor '$reactor_name' not found in $summary_file")
                continue
            end
            lcoe_mean = reactor_row[1, :mean]
        end

        push!(N_values, N_unit)
        push!(mean_lcoe, lcoe_mean)
        push!(scenario_labels, tag)
    end

    # Sort by N_values
    sort_idx = sortperm(N_values)
    N_values = N_values[sort_idx]
    mean_lcoe = mean_lcoe[sort_idx]
    scenario_labels = scenario_labels[sort_idx]

    # Create plot
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1],
              xlabel = "Number of Units Built (N)",
              ylabel = "Mean LCOE [USD/MWh]",
              title = isnothing(reactor_name) ? "Learning Curve: Mean LCOE vs Experience" : "Learning Curve: $reactor_name")

    # Add baseline as dashed reference line (SOAK level)
    baseline_file = "$outputpath/mcs-lcoe_summary-$opt_scaling-baseline.csv"
    if isfile(baseline_file)
        df_baseline = CSV.File(baseline_file) |> DataFrame
        if "mean" in names(df_baseline)
            if isnothing(reactor_name)
                # Average across all reactors
                baseline_mean = mean(filter(!ismissing, df_baseline[!, :mean]))
            else
                # Specific reactor
                reactor_row = filter(row -> row.variable == reactor_name, df_baseline)
                if nrow(reactor_row) > 0
                    baseline_mean = reactor_row[1, :mean]
                else
                    baseline_mean = nothing
                end
            end

            if !isnothing(baseline_mean)
                hlines!(ax, [baseline_mean], color = :black, linestyle = :dash, linewidth = 2.5, label = "Baseline (SOAK)")
            end
        end
    end

    # Plot learning curve line
    if !isempty(N_values)
        lines!(ax, N_values, mean_lcoe, color = :blue, linewidth = 3, label = "With Learning")
    end

    # Plot markers
    scatter!(ax, N_values, mean_lcoe, color = :blue, markersize = 15, marker = :circle)

    # Add text labels for each point
    for (i, (n, lcoe, label)) in enumerate(zip(N_values, mean_lcoe, scenario_labels))
        text!(ax, n, lcoe,
              text = "  $(round(lcoe, digits=1))",
              align = (:left, :center),
              fontsize = 12)
    end

    # Add legend
    axislegend(ax, position = :rt)

    # Add grid
    ax.xgridvisible = true
    ax.ygridvisible = true

    return fig
end

"""
    learning_curve_comparison_plot(outputpath::String, opt_scaling::String,
                                   learning_scenarios::Vector, pjs::Vector)

Create a multi-panel plot comparing learning curves for different reactor scales.

# Arguments
- `outputpath::String`: Path to directory containing CSV output files
- `opt_scaling::String`: Scaling method used
- `learning_scenarios::Vector`: Vector of tuples (N_unit, tag)
- `pjs::Vector`: Vector of project structs (to get scale information)

# Returns
- Figure object with comparison plot
"""
function learning_curve_comparison_plot(outputpath::String, opt_scaling::String,
                                       learning_scenarios::Vector, pjs::Vector)

    # Group reactors by scale
    scale_groups = Dict("Micro" => String[], "SMR" => String[], "Large" => String[])

    for pj in pjs
        if haskey(scale_groups, pj.scale)
            push!(scale_groups[pj.scale], pj.name)
        end
    end

    # Create figure with 3 subplots
    fig = Figure(size = (1800, 500))

    scale_order = ["Micro", "SMR", "Large"]
    colors = Dict("Micro" => :red, "SMR" => :blue, "Large" => :green)

    for (col_idx, scale) in enumerate(scale_order)
        if isempty(scale_groups[scale])
            continue
        end

        ax = Axis(fig[1, col_idx],
                  xlabel = "Number of Units Built (N)",
                  ylabel = "Mean LCOE [USD/MWh]",
                  title = "$scale Reactors - Learning Curve")

        # Collect data for this scale
        N_values = Int[]
        mean_lcoe = Float64[]

        for (N_unit, tag) in learning_scenarios
            summary_file = "$outputpath/mcs-lcoe_summary-$opt_scaling-$tag.csv"

            if !isfile(summary_file)
                continue
            end

            df_summary = CSV.File(summary_file) |> DataFrame

            # Calculate mean across reactors in this scale
            # In describe() output: each row is a reactor (variable column), "mean" is a column
            if !("mean" in names(df_summary))
                continue
            end

            scale_reactors = scale_groups[scale]

            # Filter to only reactors in this scale
            df_scale = filter(row -> string(row.variable) in scale_reactors, df_summary)

            if nrow(df_scale) == 0
                continue
            end

            # Get mean values for reactors in this scale
            mean_values = filter(!ismissing, df_scale[!, :mean])
            if isempty(mean_values)
                continue
            end
            lcoe_mean = mean(mean_values)

            push!(N_values, N_unit)
            push!(mean_lcoe, lcoe_mean)
        end

        # Add baseline as dashed reference line (SOAK level) for this scale
        baseline_file = "$outputpath/mcs-lcoe_summary-$opt_scaling-baseline.csv"
        if isfile(baseline_file)
            df_baseline = CSV.File(baseline_file) |> DataFrame
            if "mean" in names(df_baseline)
                scale_reactors = scale_groups[scale]
                df_scale_baseline = filter(row -> string(row.variable) in scale_reactors, df_baseline)

                if nrow(df_scale_baseline) > 0
                    baseline_mean_values = filter(!ismissing, df_scale_baseline[!, :mean])
                    if !isempty(baseline_mean_values)
                        baseline_mean = mean(baseline_mean_values)
                        hlines!(ax, [baseline_mean], color = :black, linestyle = :dash,
                               linewidth = 2, label = "Baseline (SOAK)")
                    end
                end
            end
        end

        # Sort by N
        if !isempty(N_values)
            sort_idx = sortperm(N_values)
            N_values = N_values[sort_idx]
            mean_lcoe = mean_lcoe[sort_idx]

            # Plot learning curve
            lines!(ax, N_values, mean_lcoe, color = colors[scale], linewidth = 3, label = "With Learning")
            scatter!(ax, N_values, mean_lcoe, color = colors[scale], markersize = 15)

            # Add value labels
            for (n, lcoe) in zip(N_values, mean_lcoe)
                text!(ax, n, lcoe,
                      text = "  $(round(lcoe, digits=1))",
                      align = (:left, :center),
                      fontsize = 10)
            end
        end

        # Add legend for each subplot
        axislegend(ax, position = :rt)

        ax.xgridvisible = true
        ax.ygridvisible = true
    end

    # Add overall title
    Label(fig[0, :], "Learning Curves by Reactor Scale",
          fontsize = 20, font = "Noto Sans Bold")

    return fig
end

"""
    si_plot_by_scale(si_results, title::String, pjs::Vector)

Create sensitivity index heatmap grouped by reactor scale (Micro/SMR/Large).
Three columns showing first-order (S) and total-order (ST) indices for each scale group.

# Arguments
- `si_results::DataFrame`: Sensitivity results with columns [si, var, reactor1, reactor2, ...]
- `title::String`: Plot title
- `pjs::Vector`: Vector of project objects (to get scale information)

# Returns
- Figure object with grouped heatmaps
"""
function si_plot_by_scale(si_results, title::String, pjs::Vector)

    # Group reactors by scale
    scale_groups = Dict("Micro" => String[], "SMR" => String[], "Large" => String[])

    for pj in pjs
        if haskey(scale_groups, pj.scale)
            push!(scale_groups[pj.scale], pj.name)
        end
    end

    # Filter SI results
    si_s = filter(:si => index -> index == "S", si_results)
    si_st = filter(:si => index -> index == "ST", si_results)
    xticks = si_s.var

    # Create figure with 3 scale groups side by side
    fig = Figure(size = (1800, 600))
    scale_order = ["Micro", "SMR", "Large"]
    # Use :deep colormap (same as old si_plot) for consistency across all scales
    colormap_si = :deep

    for (col_idx, scale) in enumerate(scale_order)
        reactors = scale_groups[scale]

        if isempty(reactors)
            continue
        end

        # Filter columns for this scale
        available_reactors = filter(r -> r in names(si_results), reactors)

        if isempty(available_reactors)
            continue
        end

        # Extract data for this scale
        data_s = Matrix(si_s[:, available_reactors])
        data_st = Matrix(si_st[:, available_reactors])

        yticks_scale = available_reactors

        # Create nested GridLayout for this scale group (prevents whitespace issues)
        gl = fig[1:2, col_idx] = GridLayout()

        # First-order indices (S)
        ax_s = Axis(gl[1, 1],
                    xticks = (1:length(xticks), xticks),
                    yticks = (1:length(yticks_scale), yticks_scale),
                    title = "$scale - First Order (S)")
        ax_s.xticklabelrotation = π / 3
        ax_s.xticklabelalign = (:right, :center)
        ax_s.yticklabelsize = 8  # Smaller labels for readability

        hmap_s = heatmap!(ax_s, data_s, colormap = colormap_si, colorrange = (0, 1))

        # Add text labels
        for i in 1:length(xticks), j in 1:length(yticks_scale)
            txtcolor = data_s[i, j] > 0.5 ? :white : :black
            text!(ax_s, "$(round(data_s[i,j], digits = 2))", position = (i, j),
                  color = txtcolor, fontsize = 9, align = (:center, :center))
        end

        # Total-order indices (ST)
        ax_st = Axis(gl[2, 1],
                     xticks = (1:length(xticks), xticks),
                     yticks = (1:length(yticks_scale), yticks_scale),
                     title = "$scale - Total Order (ST)")
        ax_st.xticklabelrotation = π / 3
        ax_st.xticklabelalign = (:right, :center)
        ax_st.yticklabelsize = 8  # Smaller labels for readability

        hmap_st = heatmap!(ax_st, data_st, colormap = colormap_si, colorrange = (0, 1))

        # Add text labels
        for i in 1:length(xticks), j in 1:length(yticks_scale)
            txtcolor = data_st[i, j] > 0.5 ? :white : :black
            text!(ax_st, "$(round(data_st[i,j], digits = 2))", position = (i, j),
                  color = txtcolor, fontsize = 9, align = (:center, :center))
        end

        # Add colorbar inside the GridLayout (prevents distortion)
        Colorbar(gl[1:2, 2], hmap_s, width = 12)
        colsize!(gl, 2, Auto(12))  # Keep colorbar slim
    end

    # Overall title
    Label(fig[0, :], title, fontsize = 20, font = "Noto Sans Bold", color = (:black, 0.25))

    return fig
end

"""
    shapley_plot_by_scale(shapley_results, title::String, pjs::Vector)

Create heatmap plot of Shapley sensitivity effects grouped by reactor scale.

Shapley effects differ from classical Sobol indices:
- Only ONE effect per parameter (no S/ST distinction)
- Properly accounts for parameter correlations (WACC × Construction Time)
- Effects sum to 1.0 (efficiency property)

# Arguments
- shapley_results: DataFrame with columns [var, reactor1, reactor2, ...]
- title: Plot title
- pjs: Vector of project structs for grouping by scale

# Returns
- Makie Figure with heatmaps grouped by scale (Micro/SMR/Large)
"""
function shapley_plot_by_scale(shapley_results, title::String, pjs::Vector)

    # Group reactors by scale
    scale_groups = Dict("Micro" => String[], "SMR" => String[], "Large" => String[])

    for pj in pjs
        if haskey(scale_groups, pj.scale)
            push!(scale_groups[pj.scale], pj.name)
        end
    end

    # Extract variable names (should be: wacc, construction_time, loadfactor, investment)
    xticks = shapley_results.var

    # Create figure with 3 scale groups side by side
    # Single row (unlike Sobol which has S and ST)
    fig = Figure(size = (1800, 400))
    scale_order = ["Micro", "SMR", "Large"]
    # Use :deep colormap (same as old si_plot) for consistency
    colormap_shapley = :deep

    for (col_idx, scale) in enumerate(scale_order)
        reactors = scale_groups[scale]

        if isempty(reactors)
            continue
        end

        # Filter columns for this scale
        available_reactors = filter(r -> r in names(shapley_results), reactors)

        if isempty(available_reactors)
            continue
        end

        # Extract data for this scale (only Shapley effects, no S/ST split)
        data_shapley = Matrix(shapley_results[:, available_reactors])

        yticks_scale = available_reactors

        # Create nested GridLayout for this scale group
        gl = fig[1, col_idx] = GridLayout()

        # Shapley effects heatmap
        ax_shapley = Axis(gl[1, 1],
                          xticks = (1:length(xticks), xticks),
                          yticks = (1:length(yticks_scale), yticks_scale),
                          title = "$scale - Shapley Effects")
        ax_shapley.xticklabelrotation = π / 3
        ax_shapley.xticklabelalign = (:right, :center)
        ax_shapley.yticklabelsize = 8  # Smaller labels for readability

        hmap_shapley = heatmap!(ax_shapley, data_shapley,
                                colormap = colormap_shapley,
                                colorrange = (0, 1))

        # Add text labels with values
        for i in 1:length(xticks), j in 1:length(yticks_scale)
            txtcolor = data_shapley[i, j] > 0.5 ? :white : :black
            text!(ax_shapley, "$(round(data_shapley[i,j], digits = 2))",
                  position = (i, j),
                  color = txtcolor, fontsize = 10, align = (:center, :center))
        end

        # Add colorbar inside the GridLayout
        Colorbar(gl[1, 2], hmap_shapley, width = 12)
        colsize!(gl, 2, Auto(12))  # Keep colorbar slim
    end

    # Overall title
    Label(fig[0, :], title, fontsize = 20, font = "Noto Sans Bold", color = (:black, 0.25))

    return fig
end

"""
    investment_plot_by_scale(pjs::Vector, scaling_plot::Vector)

Create investment comparison plot grouped by reactor scale (Micro/SMR/Large).
Three subplots showing manufacturer estimates vs. scaling methods for each group.

# Arguments
- `pjs::Vector`: Vector of project objects
- `scaling_plot::Vector`: Scaling parameter bounds [lower, upper]

# Returns
- Figure object with grouped investment comparisons
"""
function investment_plot_by_scale(pjs::Vector, scaling_plot::Vector)

    # Group reactors by scale
    scale_groups = Dict("Micro" => Int[], "SMR" => Int[], "Large" => Int[])

    for (idx, pj) in enumerate(pjs)
        if haskey(scale_groups, pj.scale)
            push!(scale_groups[pj.scale], idx)
        end
    end

    # Create figure with 3 subplots
    fig = Figure(size = (1800, 800))
    scale_order = ["Micro", "SMR", "Large"]
    colors = Dict("Micro" => :red, "SMR" => :blue, "Large" => :green)

    for (row_idx, scale) in enumerate(scale_order)
        indices = scale_groups[scale]

        if isempty(indices)
            continue
        end

        # Generate scaled investments for this group
        scaled_investments = DataFrame()
        reactor_names = String[]

        for idx in indices
            scaled_investments.res = gen_scaled_investment(scaling_plot, pjs[idx])
            rename!(scaled_investments, :res => pjs[idx].name)
            push!(reactor_names, pjs[idx].name)
        end

        # Create axis
        ax = Axis(fig[row_idx, 1],
                  yticks = (1:length(reactor_names), reactor_names),
                  xscale = log10,
                  xlabel = "[USD/MW]",
                  title = "$scale Reactors")

        # Plot range bars and manufacturer estimates
        rothwell = rangebars!(ax, 1.2:length(reactor_names)+0.2,
                             collect(scaled_investments[5,:]),
                             collect(scaled_investments[4,:]),
                             linewidth = 6, whiskerwidth = 12, direction = :x,
                             color = colors[scale], alpha = 0.4)

        roulstone = rangebars!(ax, 1:length(reactor_names),
                              collect(scaled_investments[3,:]),
                              collect(scaled_investments[2,:]),
                              linewidth = 6, whiskerwidth = 12, direction = :x,
                              color = colors[scale], alpha = 0.6)

        manufacturer = scatter!(ax, collect(scaled_investments[1,:]),
                               1:length(reactor_names),
                               marker = :star5, color = :black, markersize = 15)

        # Add Carelli scaling (row 6) as orange circles
        carelli = scatter!(ax, collect(scaled_investments[6,:]),
                          1:length(reactor_names),
                          marker = :circle, color = :orange, markersize = 12,
                          strokewidth = 2, strokecolor = :darkorange)

        # Add legend only to first subplot
        if row_idx == 1
            Legend(fig[row_idx, 1],
                   [manufacturer, roulstone, rothwell, carelli],
                   ["Manufacturer", "Roulstone", "Rothwell", "Carelli"],
                   tellheight = false, tellwidth = false,
                   halign = :right, valign = :top,
                   framevisible = true, orientation = :vertical)
        end
    end

    # Overall title
    Label(fig[0, :], "Investment Cost Comparison by Reactor Scale",
          fontsize = 20, font = "Noto Sans Bold")

    return fig
end

"""
    lcoe_threshold_probability_plot(lcoe_results::DataFrame, pjs::Vector;
                                    thresholds::Vector{Float64}=collect(0:10:300))

Create plot showing probability of LCOE being at or below threshold values,
grouped by reactor scale (Micro/SMR/Large).

# Arguments
- `lcoe_results::DataFrame`: LCOE Monte Carlo results (one column per reactor)
- `pjs::Vector`: Vector of project objects (to get scale information)
- `thresholds::Vector{Float64}`: LCOE threshold values to evaluate (default: 0 to 300 in steps of 10)

# Returns
- Figure object with threshold probability curves
"""
function lcoe_threshold_probability_plot(lcoe_results::DataFrame, pjs::Vector;
                                        thresholds::Vector{Float64}=collect(0:10:300))

    # Group reactors by scale
    scale_groups = Dict("Micro" => String[], "SMR" => String[], "Large" => String[])

    for pj in pjs
        if haskey(scale_groups, pj.scale) && (pj.name in names(lcoe_results))
            push!(scale_groups[pj.scale], pj.name)
        end
    end

    # Calculate probabilities for each scale
    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1, 1],
              xlabel = "LCOE Threshold [USD/MWh]",
              ylabel = "Probability LCOE > Threshold [%]",
              title = "Probability of Exceeding Cost Thresholds by Reactor Scale")

    scale_order = ["Micro", "SMR", "Large"]
    colors = Dict("Micro" => :red, "SMR" => :blue, "Large" => :green)
    markers = Dict("Micro" => :circle, "SMR" => :rect, "Large" => :diamond)

    legend_items = []
    legend_labels = String[]

    for scale in scale_order
        reactors = scale_groups[scale]

        if isempty(reactors)
            continue
        end

        # Combine all LCOE values for this scale (more efficient concatenation)
        scale_lcoe_vectors = [lcoe_results[!, reactor] for reactor in reactors]
        scale_lcoe = vcat(scale_lcoe_vectors...)

        # Calculate probabilities at each threshold (vectorized for speed)
        # P(LCOE > threshold) - probability of EXCEEDING the threshold
        n_total = length(scale_lcoe)
        probabilities = [sum(scale_lcoe .> t) / n_total * 100 for t in thresholds]

        # Plot line
        l = lines!(ax, thresholds, probabilities,
                   color = colors[scale], linewidth = 3,
                   label = "$scale (N=$(length(reactors)))")

        # Add markers
        scatter!(ax, thresholds, probabilities,
                 color = colors[scale], marker = markers[scale],
                 markersize = 8)

        push!(legend_items, l)
        push!(legend_labels, "$scale (N=$(length(reactors)))")
    end

    # Add reference lines
    hlines!(ax, [50], color = :gray, linestyle = :dash, linewidth = 2, label = "50% probability")

    # Add legend
    axislegend(ax, position = :rb)

    # Add grid
    ax.xgridvisible = true
    ax.ygridvisible = true

    return fig
end

"""
    mcs_plot_regional(lcoe_results, pjs; scale_filter="Large")

Create side-by-side percentile range plots comparing regional LCOE distributions.

# Arguments
- `lcoe_results`: DataFrame with LCOE simulation results (columns = reactor names)
- `pjs`: Vector of project objects containing reactor metadata
- `scale_filter`: Filter reactors by scale (default: "Large"). Use "All" to show all scales.

# Returns
- Figure with regional comparison using percentile ranges

Creates separate panels for each region showing:
- Light colored band: 10th-90th percentile range
- Dark colored band: 25th-75th percentile range (IQR)
- Black horizontal line: Median
- Colored diamond: Mean

Follows methodology from Weibezahn et al. (2023), OECD-NEA (2020), and
Lovering et al. (2016) using percentile bands instead of full ranges
to avoid visualization issues from extreme tails.
"""
function mcs_plot_regional(lcoe_results, pjs; scale_filter="Large")
    # Convert pjs to DataFrame for easier filtering
    pjs_dat = DataFrame(
        name = [pj.name for pj in pjs],
        scale = [pj.scale for pj in pjs],
        region = [pj.region for pj in pjs]
    )

    # Filter reactors by scale (or show all if scale_filter="All")
    if scale_filter == "All"
        scale_reactors = pjs_dat
    else
        scale_reactors = filter(row -> row.scale == scale_filter, pjs_dat)
    end

    if nrow(scale_reactors) == 0
        @warn("No reactors found with scale=$scale_filter")
        return Figure()
    end

    # Get unique regions
    regions = unique(scale_reactors.region)
    n_regions = length(regions)

    # Create figure with panels for each region
    fig = Figure(size=(500 * n_regions, 700))

    # Define colors for regions (Weibezahn et al. 2023 / OECD-NEA groupings)
    region_colors = Dict(
        "Western / Developed" => :blue,
        "Emerging Asia" => :red,
        "Eastern Europe" => :green,
        "South America" => :purple,
        "Middle East / Africa" => :orange,
        "Other" => :gray
    )

    # Plot each region
    for (col_idx, region) in enumerate(regions)
        region_reactors = filter(row -> row.region == region, scale_reactors)
        reactor_names = region_reactors.name

        if length(reactor_names) == 0
            continue
        end

        # Create axis for this region
        ax = Axis(fig[1, col_idx],
                 title="$region $(scale_filter) Reactors",
                 xlabel="Reactor",
                 ylabel="LCOE [USD/MWh]",
                 xticks=(1:length(reactor_names), reactor_names),
                 xticklabelrotation=π/4)

        ax.xticklabelsize = 10

        # Plot percentile ranges for each reactor (following Weibezahn et al. 2023, OECD-NEA)
        # Shows 10th-90th percentile (outer) and 25th-75th percentile (IQR, inner)
        for (i, reactor_name) in enumerate(reactor_names)
            if reactor_name in names(lcoe_results)
                reactor_lcoe = lcoe_results[!, reactor_name]
                color = get(region_colors, region, :gray)

                # Calculate percentiles
                p10 = quantile(reactor_lcoe, 0.10)
                p25 = quantile(reactor_lcoe, 0.25)
                p50 = quantile(reactor_lcoe, 0.50)  # median
                p75 = quantile(reactor_lcoe, 0.75)
                p90 = quantile(reactor_lcoe, 0.90)
                reactor_mean = mean(reactor_lcoe)

                # 10th-90th percentile range (outer, lighter)
                rangebars!(ax, [i], [p10], [p90],
                          color=(color, 0.3), linewidth=8,
                          whiskerwidth=8, direction=:y)

                # 25th-75th percentile range (IQR, darker)
                rangebars!(ax, [i], [p25], [p75],
                          color=(color, 0.7), linewidth=8,
                          whiskerwidth=8, direction=:y)

                # Median marker
                scatter!(ax, [i], [p50], color=:black, marker=:hline,
                        markersize=15, strokewidth=2)

                # Mean marker (diamond)
                scatter!(ax, [i], [reactor_mean], color=color, marker=:diamond,
                        markersize=12, strokecolor=:black, strokewidth=1)
            end
        end

        # Add regional mean line
        region_lcoe_vectors = [lcoe_results[!, name] for name in reactor_names if name in names(lcoe_results)]
        if !isempty(region_lcoe_vectors)
            region_all_lcoe = vcat(region_lcoe_vectors...)
            region_mean = mean(region_all_lcoe)
            hlines!(ax, [region_mean], color=:black, linestyle=:dash, linewidth=2,
                   label="Regional Mean: $(round(region_mean, digits=1)) USD/MWh")

            # Add legend
            axislegend(ax, position=:rt)
        end

        # Grid
        ax.ygridvisible = true
    end

    # Add overall title
    Label(fig[0, :], "Regional LCOE Comparison: $(scale_filter) Reactors",
          fontsize=20, font="Noto Sans Bold")

    # Add explanation subtitle
    Label(fig[0, :, Top()], "Bands show 10th-90th (light) and 25th-75th (dark) percentiles. Black line = median, diamond = mean.",
          fontsize=12, padding=(0, 0, 5, 0))

    return fig
end

"""
    mcs_plot_regional_combined(lcoe_results, pjs; scale_filter="Large")

Create combined percentile range plot comparing all regions on one axis.

# Arguments
- `lcoe_results`: DataFrame with LCOE simulation results
- `pjs`: Vector of project objects
- `scale_filter`: Filter reactors by scale (default: "Large"). Use "All" to show all scales.

# Returns
- Figure with combined regional comparison

Shows all regions on same axis for direct comparison using:
- Light colored band: 10th-90th percentile range
- Dark colored band: 25th-75th percentile range (IQR)
- Black horizontal line: Median
- Black diamond: Mean

Follows methodology from Weibezahn et al. (2023), OECD-NEA (2020), and
Lovering et al. (2016) to show meaningful ranges without extreme outliers.
"""
function mcs_plot_regional_combined(lcoe_results, pjs; scale_filter="Large")
    # Convert pjs to DataFrame
    pjs_dat = DataFrame(
        name = [pj.name for pj in pjs],
        scale = [pj.scale for pj in pjs],
        region = [pj.region for pj in pjs]
    )

    # Filter by scale (or show all if scale_filter="All")
    if scale_filter == "All"
        scale_reactors = pjs_dat
    else
        scale_reactors = filter(row -> row.scale == scale_filter, pjs_dat)
    end

    if nrow(scale_reactors) == 0
        @warn("No reactors found with scale=$scale_filter")
        return Figure()
    end

    # Get unique regions
    regions = sort(unique(scale_reactors.region))
    n_regions = length(regions)

    # Create figure
    fig = Figure(size=(1200, 700))
    ax = Axis(fig[1, 1],
             title="Regional LCOE Comparison: $(scale_filter) Reactors",
             xlabel="Region",
             ylabel="LCOE [USD/MWh]",
             xticks=(1:n_regions, regions),
             xticklabelrotation=π/6)

    # Define colors (Weibezahn et al. 2023 / OECD-NEA groupings)
    region_colors = Dict(
        "Western / Developed" => :blue,
        "Emerging Asia" => :red,
        "Eastern Europe" => :green,
        "South America" => :purple,
        "Middle East / Africa" => :orange,
        "Other" => :gray
    )

    # Collect data for each region
    regional_means = Float64[]
    regional_stds = Float64[]

    for (i, region) in enumerate(regions)
        region_reactors = filter(row -> row.region == region, scale_reactors)
        reactor_names = region_reactors.name

        # Collect all LCOE values for this region
        region_lcoe_vectors = [lcoe_results[!, name] for name in reactor_names if name in names(lcoe_results)]

        if !isempty(region_lcoe_vectors)
            region_all_lcoe = vcat(region_lcoe_vectors...)
            color = get(region_colors, region, :gray)

            # Calculate percentiles (following Weibezahn et al. 2023, OECD-NEA 2020)
            p10 = quantile(region_all_lcoe, 0.10)
            p25 = quantile(region_all_lcoe, 0.25)
            p50 = quantile(region_all_lcoe, 0.50)  # median
            p75 = quantile(region_all_lcoe, 0.75)
            p90 = quantile(region_all_lcoe, 0.90)
            region_mean = mean(region_all_lcoe)

            # 10th-90th percentile range (outer band, lighter)
            rangebars!(ax, [i], [p10], [p90],
                      color=(color, 0.3), linewidth=30,
                      whiskerwidth=15, direction=:y)

            # 25th-75th percentile range (IQR, inner band, darker)
            rangebars!(ax, [i], [p25], [p75],
                      color=(color, 0.7), linewidth=30,
                      whiskerwidth=15, direction=:y)

            # Median line
            scatter!(ax, [i], [p50], color=:black, marker=:hline,
                    markersize=25, strokewidth=3)

            # Mean marker (diamond) - will be added separately below

            # Store statistics
            push!(regional_means, region_mean)
            push!(regional_stds, std(region_all_lcoe))
        end
    end

    # Add mean markers
    if !isempty(regional_means)
        scatter!(ax, 1:length(regional_means), regional_means,
                color=:black, marker=:diamond, markersize=15,
                label="Regional Mean")

        # Add legend
        axislegend(ax, position=:rt)
    end

    # Grid
    ax.ygridvisible = true

    # Add explanation text
    Label(fig[0, :], "Bands show 10th-90th (light) and 25th-75th (dark) percentiles. Black line = median, diamond = mean.",
          fontsize=12, tellwidth=false)

    return fig
end


"""
    wacc_sensitivity_plot(outputpath, opt_scaling, pjs_dat, wacc_bin_centers)

Creates a WACC sensitivity plot showing how median LCOE varies with discount rate for different reactor scales.

**EFFICIENT VERSION:** Bins existing simulation results instead of re-running 720k simulations.

Arguments:
- `outputpath`: Path to directory containing CSV results
- `opt_scaling`: Scaling option (e.g., "roulstone")
- `pjs_dat`: DataFrame with reactor metadata (name, scale)
- `wacc_bin_centers`: Vector of WACC percentages for bin centers (e.g., 0:1:15 for 0%, 1%, ..., 15%)

Returns:
- Figure with WACC (%) vs median LCOE (USD/MWh) for Micro, SMR, and Large reactors
"""
function wacc_sensitivity_plot(outputpath, opt_scaling, pjs_dat, wacc_bin_centers)

    @info("Generating WACC sensitivity plot from existing results (0 new simulations)")

    # Read saved results
    lcoe_results = CSV.read("$outputpath/mcs-lcoe_results-$opt_scaling.csv", DataFrame)
    wacc_values = CSV.read("$outputpath/mcs-wacc_values-$opt_scaling.csv", DataFrame)

    # Storage for binned results
    wacc_percentages = Float64[]
    micro_median_lcoe = Float64[]
    smr_median_lcoe = Float64[]
    large_median_lcoe = Float64[]

    # For each WACC bin center
    for wacc_pct in wacc_bin_centers
        @info("  Binning WACC = $(wacc_pct)%")

        # Define bin edges (±0.5% around center)
        wacc_lower = (wacc_pct - 0.5) / 100.0
        wacc_upper = (wacc_pct + 0.5) / 100.0

        # Storage for this bin
        micro_lcoe_bin = Float64[]
        smr_lcoe_bin = Float64[]
        large_lcoe_bin = Float64[]

        # For each reactor
        for (idx, reactor_name) in enumerate(names(lcoe_results))
            # Get reactor scale
            scale = pjs_dat[pjs_dat.name .== reactor_name, :scale][1]

            # Get WACC and LCOE values for this reactor
            reactor_wacc = wacc_values[!, reactor_name]
            reactor_lcoe = lcoe_results[!, reactor_name]

            # Filter to this WACC bin
            in_bin = (reactor_wacc .>= wacc_lower) .& (reactor_wacc .< wacc_upper)
            lcoe_in_bin = reactor_lcoe[in_bin]

            # Append to appropriate scale bin
            if scale == "Micro"
                append!(micro_lcoe_bin, lcoe_in_bin)
            elseif scale == "SMR"
                append!(smr_lcoe_bin, lcoe_in_bin)
            elseif scale == "Large"
                append!(large_lcoe_bin, lcoe_in_bin)
            end
        end

        # Calculate median LCOE for each scale in this bin
        push!(wacc_percentages, wacc_pct)
        push!(micro_median_lcoe, isempty(micro_lcoe_bin) ? NaN : median(micro_lcoe_bin))
        push!(smr_median_lcoe, isempty(smr_lcoe_bin) ? NaN : median(smr_lcoe_bin))
        push!(large_median_lcoe, isempty(large_lcoe_bin) ? NaN : median(large_lcoe_bin))

        @info("    Bin counts: Micro=$(length(micro_lcoe_bin)), SMR=$(length(smr_lcoe_bin)), Large=$(length(large_lcoe_bin))")
    end

    # Create plot
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1,1],
             xlabel="Discount rate (WACC, %)",
             ylabel="Median LCOE [USD/MWh]",
             title="LCOE Sensitivity to Discount Rate by Reactor Scale")

    # Plot lines for each scale
    micro_line = lines!(ax, wacc_percentages, micro_median_lcoe,
                       label="Micro", linewidth=3, color=:blue)
    smr_line = lines!(ax, wacc_percentages, smr_median_lcoe,
                     label="SMR", linewidth=3, color=:green)
    large_line = lines!(ax, wacc_percentages, large_median_lcoe,
                       label="Large", linewidth=3, color=:red)

    # Add markers
    scatter!(ax, wacc_percentages, micro_median_lcoe,
            color=:blue, markersize=10)
    scatter!(ax, wacc_percentages, smr_median_lcoe,
            color=:green, markersize=10)
    scatter!(ax, wacc_percentages, large_median_lcoe,
            color=:red, markersize=10)

    # Add legend
    Legend(fig[1,2], ax, framevisible=true)

    # Grid
    ax.xgridvisible = true
    ax.ygridvisible = true

    @info("WACC sensitivity plot complete (used existing data, 0 new simulations)")

    return fig
end

"""
    learning_curve_plot(data_path::String, opt_scaling::String)

Create learning curve plots showing LCOE vs cumulative units (N) for different learning rates.

# Arguments
- `data_path`: Path to directory containing learning curve CSV data
- `opt_scaling`: Scaling method used (e.g., "rothwell", "roulstone")

# Returns
- Tuple of 3 Figure objects: (Micro, SMR, Large) learning curve plots

Each plot shows:
- X-axis: Cumulative units (N)
- Y-axis: LCOE [USD/MWh]
- Lines: One per reactor × learning rate combination
- Shaded bands: P10-P90 uncertainty ranges
- Colors: By reactor type (PWR/BWR=orangered, HTR=gold, SFR=teal)
"""
function learning_curve_plot(data_path::String, opt_scaling::String)
    
    @info("Generating learning curve plots")
    
    # Load learning results
    csv_file = "$data_path/learning-lcoe_curves-$opt_scaling.csv"
    if !isfile(csv_file)
        @warn("Learning curve data not found: $csv_file")
        @warn("Run run_4_learning.jl first to generate learning curve data")
        return (nothing, nothing, nothing)
    end
    
    df = CSV.read(csv_file, DataFrame)
    
    # Define scales and learning rate labels
    scales = ["Micro", "SMR", "Large"]
    lr_labels = Dict(0.05 => "Pessimistic (5%)", 0.10 => "Base (10%)", 0.15 => "Optimistic (15%)")
    lr_colors = Dict(0.05 => :red, 0.10 => :blue, 0.15 => :green)
    
    # Color scheme by reactor type
    type_colors = Dict("PWR" => :orangered, "BWR" => :orangered, "HTR" => :gold, "SFR" => :teal)
    
    figs = []
    
    for scale in scales
        @info("  Creating learning curve for $scale reactors")
        
        df_scale = filter(row -> row.scale == scale, df)
        
        if nrow(df_scale) == 0
            @warn("No data for $scale scale, skipping")
            push!(figs, nothing)
            continue
        end
        
        fig = Figure(size=(1400, 900))
        
        # Get unique reactors for this scale
        reactors = unique(df_scale.reactor)
        
        # Determine grid layout based on number of reactors
        n_reactors = length(reactors)
        n_cols = min(3, n_reactors)  # Max 3 columns
        n_rows = ceil(Int, n_reactors / n_cols)
        
        for (reactor_idx, reactor) in enumerate(reactors)
            df_reactor = filter(row -> row.reactor == reactor, df_scale)
            
            # Calculate subplot position
            row_idx = div(reactor_idx - 1, n_cols) + 1
            col_idx = mod(reactor_idx - 1, n_cols) + 1
            
            # Create axis for this reactor
            ax = Axis(fig[row_idx, col_idx],
                xlabel = "Cumulative Units (N)",
                ylabel = "LCOE [USD/MWh]",
                title = reactor)
            
            # Get reactor type for coloring
            reactor_type = first(df_reactor.type)
            base_color = get(type_colors, reactor_type, :gray)
            
            # Plot each learning rate
            for (lr_idx, LR) in enumerate([0.05, 0.10, 0.15])
                df_lr = filter(row -> row.LR == LR, df_reactor)
                
                if nrow(df_lr) == 0
                    continue
                end
                
                # Sort by N_unit
                sort!(df_lr, :N_unit)
                
                # Adjust line style and opacity by learning rate
                line_alpha = lr_idx == 2 ? 1.0 : 0.7  # Base case (10%) most prominent
                line_width = lr_idx == 2 ? 3 : 2
                
                # Plot median line
                lines!(ax, df_lr.N_unit, df_lr.LCOE_median,
                    label = lr_labels[LR],
                    color = (base_color, line_alpha),
                    linewidth = line_width)
                
                # Add uncertainty band (P10-P90)
                band!(ax, df_lr.N_unit, df_lr.LCOE_p10, df_lr.LCOE_p90,
                    color = (base_color, 0.15))
                
                # Add markers at key points
                scatter!(ax, df_lr.N_unit, df_lr.LCOE_median,
                    color = (base_color, line_alpha),
                    markersize = 6)
            end
            
            # Add legend to first subplot only
            if reactor_idx == 1
                axislegend(ax, position=:rt, framevisible=true)
            end
            
            # Enable grid
            ax.xgridvisible = true
            ax.ygridvisible = true
        end
        
        # Add overall title
        Label(fig[0, :], "Learning Curves: $scale Reactors ($opt_scaling scaling)",
              fontsize = 20, font = "Noto Sans Bold")
        
        push!(figs, fig)
    end
    
    @info("Learning curve plots generated")
    
    return (figs[1], figs[2], figs[3])  # Return as tuple (Micro, SMR, Large)
end

"""
    lcoe_threshold_probability_plot_styled(lcoe_results::DataFrame, pjs::Vector;
                                          thresholds::Vector{Float64}=collect(0:10:300))

Create a styled LCOE threshold probability plot showing individual reactor curves.

# Arguments
- `lcoe_results`: DataFrame with LCOE simulation results (columns = reactor names)
- `pjs`: Vector of project objects containing reactor metadata
- `thresholds`: Vector of LCOE thresholds to evaluate (default: 0-300 USD/MWh)

# Returns
- Figure with individual reactor probability curves colored by scale

Creates a single panel showing:
- All individual reactor probability curves
- Color palettes by scale (with gradients for visual separation):
  * Micro: orange/coral tones (#FF6B6B to #FF8E53)
  * SMR: blue/cyan tones (#4ECDC4 to #44A8D8)
  * Large: green/cyan tones (#95E1D3 to #48C9B0)
- 50% probability reference line (dashed gray)
- Grouped legend by scale with reactor counts
- Transparent lines (alpha=0.7) for better overlapping visibility
- High resolution (1400x800) for thesis quality
"""
function lcoe_threshold_probability_plot_styled(lcoe_results::DataFrame, pjs::Vector;
                                               thresholds::Vector{Float64}=collect(0:10:300))
    
    # Define color palettes by scale (RGB tuples normalized to 0-1)
    # Converted from hex: #RRGGBB -> (R/255, G/255, B/255)
    scale_palettes = Dict(
        "Micro" => [
            (1.0, 0.42, 0.42),  # #FF6B6B - Coral red
            (1.0, 0.48, 0.48),  # #FF7B7B
            (1.0, 0.55, 0.42),  # #FF8B6B
            (1.0, 0.56, 0.33)   # #FF8E53 - Orange
        ],
        "SMR" => [
            (0.31, 0.80, 0.77),  # #4ECDC4 - Cyan
            (0.29, 0.74, 0.78),  # #4BBDC8
            (0.28, 0.68, 0.82),  # #48AED0
            (0.27, 0.66, 0.85)   # #44A8D8 - Blue
        ],
        "Large" => [
            (0.58, 0.88, 0.83),  # #95E1D3 - Light cyan
            (0.50, 0.86, 0.78),  # #7FDCC8
            (0.40, 0.83, 0.74),  # #66D3BC
            (0.28, 0.79, 0.69)   # #48C9B0 - Teal green
        ]
    )
    
    # Group reactors by scale
    scale_groups = Dict("Micro" => String[], "SMR" => String[], "Large" => String[])
    
    for pj in pjs
        if haskey(scale_groups, pj.scale) && (pj.name in names(lcoe_results))
            push!(scale_groups[pj.scale], pj.name)
        end
    end
    
    # Create figure with thesis-quality resolution
    fig = Figure(size = (1400, 800))
    
    ax = Axis(fig[1, 1],
              xlabel = "LCOE Threshold [USD/MWh]",
              ylabel = "Probability of Exceeding Threshold [%]",
              xlabelsize = 16,
              ylabelsize = 16,
              xticklabelsize = 14,
              yticklabelsize = 14)
    
    # Enable grid with subtle styling
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridcolor = (:gray, 0.3)
    ax.ygridcolor = (:gray, 0.3)
    ax.xgridwidth = 0.5
    ax.ygridwidth = 0.5
    
    scale_order = ["Micro", "SMR", "Large"]
    legend_elements = []
    legend_labels = String[]
    
    # Plot individual reactor curves
    for scale in scale_order
        reactors = scale_groups[scale]

        if isempty(reactors)
            continue
        end

        n_reactors = length(reactors)
        palette = scale_palettes[scale]

        # Generate colors for this scale (interpolate if needed)
        if n_reactors <= length(palette)
            colors = palette[1:n_reactors]
        else
            # Interpolate colors if more reactors than palette colors
            colors = [palette[max(1, min(length(palette), floor(Int, (i-1)/(n_reactors-1) * (length(palette)-1)) + 1))]
                     for i in 1:n_reactors]
        end

        # Plot each reactor and add to legend
        for (idx, reactor) in enumerate(reactors)
            lcoe_data = lcoe_results[!, reactor]

            # Calculate probabilities: P(LCOE > threshold)
            n_total = length(lcoe_data)
            probabilities = [sum(lcoe_data .> t) / n_total * 100 for t in thresholds]

            # Plot line with transparency and save for legend
            line = lines!(ax, thresholds, probabilities,
                         color = (colors[idx], 0.7),
                         linewidth = 2.5)

            # Add to legend
            push!(legend_elements, line)
            push!(legend_labels, reactor)
        end
    end

    # Add 50% probability reference line
    hlines!(ax, [50.0],
            color = :gray,
            linestyle = :dash,
            linewidth = 2,
            label = "50% threshold")

    # Add legend with individual reactor names (2 columns for better space usage)
    Legend(fig[1, 2],
           legend_elements,
           legend_labels,
           framevisible = true,
           labelsize = 11,
           nbanks = 2,  # Two columns for 25 reactors
           titlesize = 13,
           patchsize = (25, 15))  # Adjust patch size for better visibility
    
    return fig
end
