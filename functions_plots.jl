"""
The investment_plot function generates a comparison plot of investment values for manufacturer vs. estimated costs for multiple projects. The function takes in two arguments: a vector pjs containing project objects to be compared, and a boolean scaling_plot which indicates whether to scale the investment values by the project capacity.
The function first generates scaled investment values for all projects using the gen_scaled_investment function, and stores them in a DataFrame. The DataFrame columns are named after the project names.
The function then creates a new figure object and defines the x-axis label and y-axis tick labels based on the input DataFrame.
The function then plots the investment values using three different plot types: scatter for the manufacturer's investment values, and rangebars for the other two project's investment values. The rangebars represent the minimum and maximum investment values for each project.
Finally, the function adds a legend to the plot and returns the figure object.
"""
function investment_plot(pjs, scaling_plot)

    scaled_investments = DataFrame();

    # generate bounds for scaled investment values for all projects
    for p in eachindex(pjs)
        scaled_investments.res = gen_scaled_investment(scaling_plot, pjs[p])
        rename!(scaled_investments,:res => pjs[p].name)
    end

    # define plot
    xlabel = "[USD/MW]";
    yticks = names(scaled_investments);

    fig_invest_comparison = Figure();

    ax_invest = Axis(fig_invest_comparison[1,1], yticks = (1:length(yticks), yticks), xscale = log10, xlabel = xlabel);

    rothwell = rangebars!(ax_invest, 1.2:length(yticks)+0.2, collect(scaled_investments[5,:]), collect(scaled_investments[4,:]), linewidth = 6, whiskerwidth = 12, direction = :x, transparency = :true, color = :green)
    roulstone = rangebars!(ax_invest, 1:length(yticks), collect(scaled_investments[3,:]), collect(scaled_investments[2,:]), linewidth = 6, whiskerwidth = 12, direction = :x, transparency = :true, color = :darkblue)
    manufacturer = scatter!(ax_invest, collect(scaled_investments[1,:]), 1:length(yticks), marker = :star5, color = :red)

    Legend(fig_invest_comparison[1, 1],
        [manufacturer, roulstone, rothwell],
        ["Manufacturer", "Roulstone", "Rothwell"],
        tellheight = false,
        tellwidth = false,
        halign = :right, valign = :bottom,
        framevisible = false, orientation = :vertical)

    return fig_invest_comparison

end

"""
The hist_invest_plot function generates a histogram plot of the investment values of a given project for a specified number of trials. The function takes in six arguments:
- `n::Int64`: The number of trials.
- `wacc::Vector`: A vector of weighted average cost of capital values for each trial (not used for the plot but needed for the gen_rand_vars function).
- `electricity_price::Vector`: A vector of electricity prices for each trial (not used for the plot but needed for the gen_rand_vars function).
- `pj::project`: An instance of a project class that contains information about the project.
- `i::Int64`: The row index of the subplot where the histogram plot will be added.
- `j::Int64`: The column index of the subplot where the histogram plot will be added.
- `hist_invest = Figure()`: The histogram plot is added to this figure object. If not provided, a new figure object will be created.

The function first generates random variables using the gen_rand_vars function for two different options of scaling. It then adds the histograms to the plot. Each histogram is normalized to show the probability density function.
"""
function hist_invest_plot(n::Int64, wacc::Vector, electricity_price::Vector, pj::project, i::Int64, j::Int64, hist_invest = Figure())

    random_vars_roulstone = gen_rand_vars(opts_scaling[2], n, wacc, electricity_price, pj)
    random_vars_rothwell = gen_rand_vars(opts_scaling[3], n, wacc, electricity_price, pj)
    # random_vars_uniform = gen_rand_vars(opts_scaling[4], n, wacc, electricity_price, pj)

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
- `electricity_price::Vector`: A vector of electricity prices for each trial (not used for the plot but needed for the gen_rand_vars function).
- `pj::project`: An instance of a project class that contains information about the project.
- `i::Int64`: The row index of the subplot where the histogram plot will be added.
- `j::Int64`: The column index of the subplot where the histogram plot will be added.
- `density_invest = Figure()`: The histogram plot is added to this figure object. If not provided, a new figure object will be created.

The function first generates random variables using the gen_rand_vars function for two different options of scaling. It then adds the probability density function to the plot.
"""
function density_invest_plot(n::Int64, wacc::Vector, electricity_price::Vector, pj::project, i::Int64, j::Int64, density_invest = Figure())

    random_vars_roulstone = gen_rand_vars(opts_scaling[2], n, wacc, electricity_price, pj)
    random_vars_rothwell = gen_rand_vars(opts_scaling[3], n, wacc, electricity_price, pj)
    # random_vars_uniform = gen_rand_vars(opts_scaling[4], n, wacc, electricity_price, pj)

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
The function generates a box plot for each of the three categories of reactor types. The x-axis of each box plot is labeled with the parameter names, and the y-axis is labeled with the ylabel argument. The color of the box plots depends on the reactor type. The figure has a title specified by the title argument, and each box plot has a label specified by the corresponding reactor type.
"""
function mcs_plot(mcs_results, title::String, ylabel::String)

    xticks_wc = names(mcs_results)[1:9];
    xticks_ht = names(mcs_results)[10:12];
    xticks_sf = names(mcs_results)[13:15];

    mcs_boxplot = Figure();

    ax_wc = Axis(mcs_boxplot[1,1], xticks = (1:length(xticks_wc), xticks_wc), ylabel = ylabel);
    ax_wc.xticklabelrotation = π / 3;
    ax_wc.yticklabelrotation = π / 2;
    ax_wc.xticklabelalign = (:right, :center);

    for i in 1:length(xticks_wc)
        boxplot!(ax_wc, fill(i,n), mcs_results[!,i], color = :green)
    end;

    Label(mcs_boxplot[1, 1, Top()], "BWR & PWR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0))

    ax_ht = Axis(mcs_boxplot[1,2], xticks = (1:length(xticks_ht), xticks_ht));
    ax_ht.xticklabelrotation = π / 3;
    ax_ht.yticklabelrotation = π / 2;
    ax_ht.xticklabelalign = (:right, :center);

    for i in 1:length(xticks_ht)
        boxplot!(ax_ht, fill(i,n), mcs_results[!,i+length(xticks_wc)], color = :red)
    end;

    Label(mcs_boxplot[1, 2, Top()], "HTR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    ax_sf = Axis(mcs_boxplot[1,3], xticks = (1:length(xticks_sf), xticks_sf));
    ax_sf.xticklabelrotation = π / 3;
    ax_sf.yticklabelrotation = π / 2;
    ax_sf.xticklabelalign = (:right, :center);

    for i in 1:length(xticks_sf)
        boxplot!(ax_sf, fill(i,n), mcs_results[!,i+length(xticks_wc)+length(xticks_sf)], color = :orange)
    end;

    Label(mcs_boxplot[1, 3, Top()], "SFR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    Label(mcs_boxplot[0, :], title, fontsize = 24, font = "Noto Sans Bold", color = (:black, 0.25))

    colsize!(mcs_boxplot.layout, 1, Relative(6/10));
    colsize!(mcs_boxplot.layout, 2, Relative(2/10));
    colsize!(mcs_boxplot.layout, 3, Relative(2/10));

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
    fig = Figure(resolution = (1800, 500))

    # Define colors for each scale
    colors = Dict("Micro" => :red, "SMR" => :blue, "Large" => :green)

    # Define scale order and column positions
    scale_order = ["Micro", "SMR", "Large"]

    for (col_idx, scale) in enumerate(scale_order)
        if haskey(lcoe_by_scale, scale) && !isempty(lcoe_by_scale[scale])
            data = lcoe_by_scale[scale]

            # Create axis for this subplot
            ax = Axis(fig[1, col_idx],
                      xlabel = "LCOE [USD/MWh]",
                      ylabel = "Probability Density",
                      title = "$scale Reactors")

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
            mean_col = df_summary[df_summary.variable .== "mean", :]
            # Get all numeric columns (skip 'variable' column)
            numeric_cols = names(mean_col)[2:end]

            # Apply scale filter if specified
            if !isnothing(scale_filter)
                # Read reactor data to get scale information
                # This is a simplified approach - you may need to adjust based on your data structure
                numeric_cols = filter(col -> occursin(scale_filter, col), numeric_cols)
            end

            mean_values = [mean_col[1, col] for col in numeric_cols if !ismissing(mean_col[1, col])]
            lcoe_mean = mean(mean_values)
        else
            # Specific reactor
            if !(reactor_name in names(df_summary))
                @warn("Reactor '$reactor_name' not found in $summary_file")
                continue
            end
            mean_row = df_summary[df_summary.variable .== "mean", :]
            lcoe_mean = mean_row[1, reactor_name]
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
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1],
              xlabel = "Number of Units Built (N)",
              ylabel = "Mean LCOE [USD/MWh]",
              title = isnothing(reactor_name) ? "Learning Curve: Mean LCOE vs Experience" : "Learning Curve: $reactor_name")

    # Plot line
    lines!(ax, N_values, mean_lcoe, color = :blue, linewidth = 3, label = "Mean LCOE")

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
    fig = Figure(resolution = (1800, 500))

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
            mean_row = df_summary[df_summary.variable .== "mean", :]
            scale_reactors = scale_groups[scale]

            # Filter to only reactors in this scale that exist in the data
            available_reactors = filter(r -> r in names(mean_row), scale_reactors)

            if isempty(available_reactors)
                continue
            end

            mean_values = [mean_row[1, col] for col in available_reactors if !ismissing(mean_row[1, col])]
            lcoe_mean = mean(mean_values)

            push!(N_values, N_unit)
            push!(mean_lcoe, lcoe_mean)
        end

        # Sort by N
        if !isempty(N_values)
            sort_idx = sortperm(N_values)
            N_values = N_values[sort_idx]
            mean_lcoe = mean_lcoe[sort_idx]

            # Plot
            lines!(ax, N_values, mean_lcoe, color = colors[scale], linewidth = 3)
            scatter!(ax, N_values, mean_lcoe, color = colors[scale], markersize = 15)

            # Add value labels
            for (n, lcoe) in zip(N_values, mean_lcoe)
                text!(ax, n, lcoe,
                      text = "  $(round(lcoe, digits=1))",
                      align = (:left, :center),
                      fontsize = 10)
            end
        end

        ax.xgridvisible = true
        ax.ygridvisible = true
    end

    # Add overall title
    Label(fig[0, :], "Learning Curves by Reactor Scale",
          fontsize = 20, font = "Noto Sans Bold")

    return fig
end
