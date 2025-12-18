"""
    decompose_capital_to_occ_idc(capital_pv_values, wacc_values, construction_time_values)

Illustrative decomposition of capital cost (from DCF) into OCC and IDC components
using the Wealer explicit IDC formula.

This is for visualization purposes only. The actual simulation uses DCF framework
where IDC is implicitly captured through present value discounting. This function
back-calculates what OCC and IDC would be if using the explicit IDC approach.

# Arguments
- `capital_pv_values`: Array of capital cost (nominal) from DCF simulation
- `wacc_values`: Array of WACC values used in simulation
- `construction_time_values`: Array of construction times used in simulation

# Returns
- `(occ_values, idc_values)`: Illustrative split of capital cost

# Formula
IDC = OCC × [(r/2)×T + (r²/6)×T²]  (Wealer formula)
TCC = OCC + IDC
Therefore: OCC = capital_pv / (1 + IDC_factor)
"""
function decompose_capital_to_occ_idc(capital_pv_values, wacc_values, construction_time_values)
    n = length(capital_pv_values)
    occ_illustrative = zeros(Float64, n)
    idc_illustrative = zeros(Float64, n)

    for i in 1:n
        # Calculate IDC factor using Wealer formula
        T = construction_time_values[i]
        r = wacc_values[i]
        idc_factor = (r/2) * T + (r^2/6) * T^2

        # capital_pv is the nominal total capital cost (TCC) in our DCF framework
        # TCC = OCC × (1 + idc_factor)
        # Solving for OCC:
        occ_illustrative[i] = capital_pv_values[i] / (1 + idc_factor)
        idc_illustrative[i] = capital_pv_values[i] - occ_illustrative[i]
    end

    return (occ_illustrative, idc_illustrative)
end

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
#     xlabel = "[EUR2025/MW]";
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

    # Increased height to accommodate rotated labels
    mcs_boxplot = Figure(size=(1600, 1100))
    n = nrow(mcs_results)  # Number of Monte Carlo samples

    # Filter data to 5-95% percentile range to remove extreme outliers
    # This applies the filtering to the DATA, not just the axis limits
    filtered_results = copy(mcs_results)

    for col in names(filtered_results)
        col_data = filtered_results[!, col]
        q05 = quantile(col_data, 0.05)
        q95 = quantile(col_data, 0.95)

        # Filter: keep only values within 5-95% percentile
        mask = (col_data .>= q05) .& (col_data .<= q95)
        filtered_results[!, col] = ifelse.(mask, col_data, missing)
    end

    # Remove rows that are all missing (if any)
    filtered_results = dropmissing(filtered_results, names(filtered_results))

    n_filtered = nrow(filtered_results)
    @info "Box plot data filtered: $n original samples → $n_filtered samples (5-95% percentile)"

    ax_wc = Axis(mcs_boxplot[1,1],
                 xticks = (1:length(xticks_wc), xticks_wc),
                 ylabel = ylabel,
                 xticklabelspace = 120.0);  # Increase space for rotated labels
    ax_wc.xticklabelrotation = π / 3;
    ax_wc.yticklabelrotation = π / 2;
    ax_wc.xticklabelalign = (:right, :center);
    ax_wc.ygridvisible = true;  # Add grid for easier reading

    for (i, reactor_name) in enumerate(xticks_wc)
        if hasproperty(filtered_results, reactor_name)
            data_vec = collect(skipmissing(filtered_results[!, reactor_name]))
            boxplot!(ax_wc, fill(i, length(data_vec)), data_vec, color = :orangered, show_outliers=false)
        else
            @warn "Reactor $reactor_name not found in MCS results"
        end
    end;

    Label(mcs_boxplot[1, 1, Top()], "BWR & PWR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0))

    ax_ht = Axis(mcs_boxplot[1,2],
                 xticks = (1:length(xticks_ht), xticks_ht),
                 xticklabelspace = 120.0);  # Increase space for rotated labels
    ax_ht.xticklabelrotation = π / 3;
    ax_ht.yticklabelrotation = π / 2;
    ax_ht.xticklabelalign = (:right, :center);
    ax_ht.ygridvisible = true;

    for (i, reactor_name) in enumerate(xticks_ht)
        if hasproperty(filtered_results, reactor_name)
            data_vec = collect(skipmissing(filtered_results[!, reactor_name]))
            boxplot!(ax_ht, fill(i, length(data_vec)), data_vec, color = :gold, show_outliers=false)
        else
            @warn "Reactor $reactor_name not found in MCS results"
        end
    end;

    Label(mcs_boxplot[1, 2, Top()], "HTR types", font = "Noto Sans Bold", padding = (0, 6, 6, 0));

    ax_sf = Axis(mcs_boxplot[1,3],
                 xticks = (1:length(xticks_sf), xticks_sf),
                 xticklabelspace = 120.0);  # Increase space for rotated labels
    ax_sf.xticklabelrotation = π / 3;
    ax_sf.yticklabelrotation = π / 2;
    ax_sf.xticklabelalign = (:right, :center);
    ax_sf.ygridvisible = true;

    for (i, reactor_name) in enumerate(xticks_sf)
        if hasproperty(filtered_results, reactor_name)
            data_vec = collect(skipmissing(filtered_results[!, reactor_name]))
            boxplot!(ax_sf, fill(i, length(data_vec)), data_vec, color = :teal, show_outliers=false)
        else
            @warn "Reactor $reactor_name not found in MCS results"
        end
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

    # Color scheme matching learning rate plots (by reactor type)
    type_colors = Dict(
        "PWR" => RGBf(0.8, 0.2, 0.2),    # Orangered
        "BWR" => RGBf(0.9, 0.5, 0.0),    # Orange
        "HTR" => RGBf(0.85, 0.7, 0.0),   # Gold
        "SFR" => RGBf(0.0, 0.5, 0.5),    # Teal
        "MSR" => RGBf(0.6, 0.0, 0.6)     # Purple
    )

    # Create figure with 3 subplots (1 row, 3 columns)
    fig = Figure(size = (2000, 700), backgroundcolor = :white)

    # Define scale order and column positions
    scale_order = ["Micro", "SMR", "Large"]

    for (col_idx, scale) in enumerate(scale_order)
        # Get all reactors for this scale
        scale_reactors = filter(pj -> pj.scale == scale, pjs)

        if isempty(scale_reactors)
            continue
        end

        # Aggregate all data for this scale
        all_scale_data = Float64[]
        for pj in scale_reactors
            if pj.name in names(lcoe_results)
                append!(all_scale_data, lcoe_results[!, pj.name])
            end
        end

        if isempty(all_scale_data)
            continue
        end

        # Filter to 5-95% percentile to remove extreme outliers
        q05 = quantile(all_scale_data, 0.05)
        q95 = quantile(all_scale_data, 0.95)
        filtered_data = filter(x -> q05 <= x <= q95, all_scale_data)

        # Create axis
        ax = Axis(fig[1, col_idx],
                  xlabel = "LCOE [EUR2025/MWh]",
                  ylabel = col_idx == 1 ? "Probability Density" : "",
                  title = "$scale Reactors",
                  titlesize = 18,
                  xlabelsize = 14,
                  ylabelsize = 14)

        # Plot separate histogram for each reactor type within this scale
        # to show the distribution while maintaining visual clarity
        type_data = Dict{String, Vector{Float64}}()

        for pj in scale_reactors
            if pj.name in names(lcoe_results)
                reactor_lcoe = filter(x -> q05 <= x <= q95, lcoe_results[!, pj.name])
                if !haskey(type_data, pj.type)
                    type_data[pj.type] = Float64[]
                end
                append!(type_data[pj.type], reactor_lcoe)
            end
        end

        # Plot stacked/overlaid histograms by type
        for (reactor_type, data) in sort(collect(type_data))
            if !isempty(data)
                color = get(type_colors, reactor_type, RGBf(0.5, 0.5, 0.5))
                hist!(ax, data,
                      bins = 40,
                      normalization = :pdf,
                      color = (color, 0.5),
                      strokewidth = 0.5,
                      strokecolor = (color, 0.8),
                      label = reactor_type)
            end
        end

        # Calculate and plot mean and median for the entire scale
        mean_val = mean(filtered_data)
        median_val = median(filtered_data)

        vlines!(ax, [mean_val],
                color = :black,
                linestyle = :solid,
                linewidth = 2.5,
                label = "Mean: $(round(mean_val, digits=1))")

        vlines!(ax, [median_val],
                color = (:black, 0.7),
                linestyle = :dash,
                linewidth = 2.5,
                label = "Median: $(round(median_val, digits=1))")

        # Add grid for readability
        ax.xgridvisible = true
        ax.ygridvisible = true
        ax.xgridcolor = (:black, 0.1)
        ax.ygridcolor = (:black, 0.1)

        # Add legend - always use right-top to avoid overlapping with data
        axislegend(ax, position = :rt, framevisible = true,
                   labelsize = 11, bgcolor = (:white, 0.9))
    end

    # Add overall title
    Label(fig[0, :], "LCOE Distribution by Reactor Scale (5-95% Percentile)",
          fontsize = 22, font = "Noto Sans Bold", color = (:black, 0.7))

    return fig
end

"""
    lcoe_scale_histogram_single(lcoe_results::DataFrame, pjs::Vector, scale::String)

Create a single histogram showing LCOE distribution for one reactor scale.
Returns a standalone figure suitable for individual PDF export.

# Arguments
- `lcoe_results::DataFrame`: DataFrame containing LCOE results from Monte Carlo simulation
- `pjs::Vector`: Vector of project objects containing reactor information
- `scale::String`: Scale to plot ("Micro", "SMR", or "Large")

# Returns
- Figure object with single histogram
"""
function lcoe_scale_histogram_single(lcoe_results::DataFrame, pjs::Vector, scale::String)

    # Color scheme matching learning rate plots (by reactor type)
    type_colors = Dict(
        "PWR" => RGBf(0.8, 0.2, 0.2),
        "BWR" => RGBf(0.9, 0.5, 0.0),
        "HTR" => RGBf(0.85, 0.7, 0.0),
        "SFR" => RGBf(0.0, 0.5, 0.5),
        "MSR" => RGBf(0.6, 0.0, 0.6)
    )

    # Create single figure - wider for better readability
    fig = Figure(size = (800, 500), backgroundcolor = :white)

    # Get all reactors for this scale
    scale_reactors = filter(pj -> pj.scale == scale, pjs)

    if isempty(scale_reactors)
        @warn "No reactors found for scale: $scale"
        return fig
    end

    # Aggregate all data for this scale
    all_scale_data = Float64[]
    for pj in scale_reactors
        if pj.name in names(lcoe_results)
            append!(all_scale_data, lcoe_results[!, pj.name])
        end
    end

    if isempty(all_scale_data)
        @warn "No LCOE data found for scale: $scale"
        return fig
    end

    # Filter to 5-95% percentile to remove extreme outliers
    q05 = quantile(all_scale_data, 0.05)
    q95 = quantile(all_scale_data, 0.95)
    filtered_data = filter(x -> q05 <= x <= q95, all_scale_data)

    # Create axis
    ax = Axis(fig[1, 1],
              xlabel = "LCOE [EUR₂₀₂₅/MWh]",
              ylabel = "Probability Density",
              title = "$scale Reactors: LCOE Distribution",
              titlesize = 18,
              xlabelsize = 14,
              ylabelsize = 14)

    # Plot separate histogram for each reactor type within this scale
    type_data = Dict{String, Vector{Float64}}()

    for pj in scale_reactors
        if pj.name in names(lcoe_results)
            reactor_lcoe = filter(x -> q05 <= x <= q95, lcoe_results[!, pj.name])
            if !haskey(type_data, pj.type)
                type_data[pj.type] = Float64[]
            end
            append!(type_data[pj.type], reactor_lcoe)
        end
    end

    # Plot stacked/overlaid histograms by type
    for (reactor_type, data) in sort(collect(type_data))
        if !isempty(data)
            color = get(type_colors, reactor_type, RGBf(0.5, 0.5, 0.5))
            hist!(ax, data,
                  bins = 40,
                  normalization = :pdf,
                  color = (color, 0.5),
                  strokewidth = 0.5,
                  strokecolor = (color, 0.8),
                  label = reactor_type)
        end
    end

    # Calculate and plot mean and median for the entire scale
    mean_val = mean(filtered_data)
    median_val = median(filtered_data)

    vlines!(ax, [mean_val],
            color = :black,
            linestyle = :solid,
            linewidth = 2.5,
            label = "Mean: $(round(mean_val, digits=1))")

    vlines!(ax, [median_val],
            color = (:black, 0.7),
            linestyle = :dash,
            linewidth = 2.5,
            label = "Median: $(round(median_val, digits=1))")

    # Add grid for readability
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridcolor = (:black, 0.1)
    ax.ygridcolor = (:black, 0.1)

    # Add legend - always use right-top to avoid overlapping with data
    legend_pos = :rt
    axislegend(ax, position = legend_pos, framevisible = true,
               labelsize = 11, bgcolor = (:white, 0.9))

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
              ylabel = "Mean LCOE [EUR2025/MWh]",
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
                  ylabel = "Mean LCOE [EUR2025/MWh]",
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

        # Add legend for each subplot (only if there are labeled elements)
        try
            axislegend(ax, position = :rt)
        catch e
            # Skip legend if no labeled plot elements exist
            @debug "No legend items for scale $scale: $e"
        end

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

    # Map variable names to display names (matching Shapley plot)
    var_display_names = Dict(
        "wacc" => "WACC",
        "construction_time" => "Construction Time",
        "capacity factor" => "Capacity Factor",
        "scaling" => "Investment OCC"
    )

    xticks_raw = si_s.var
    xticks = [get(var_display_names, var, var) for var in xticks_raw]

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

    # Extract variable names and map to display names
    var_display_names = Dict(
        "wacc" => "WACC",
        "construction_time" => "Construction Time",
        "capacity factor" => "Capacity Factor",
        "scaling" => "Investment OCC"
    )

    xticks_raw = shapley_results.var
    xticks = [get(var_display_names, var, var) for var in xticks_raw]

    # Create figure with 3 scale groups side by side
    # Single row (unlike Sobol which has S and ST)
    # Increased width to 2200 to give more space to SMR panel
    # Increased height to 650 for better reactor spacing and readability
    fig = Figure(size = (2200, 650))
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

        # Identify reactor types for color coding
        reactor_types_in_scale = []
        type_indicator_colors_map = Dict(
            "PWR" => RGBf(0.9, 0.5, 0.1),    # Orange
            "BWR" => RGBf(0.9, 0.5, 0.1),    # Orange (same as PWR)
            "HTR" => RGBf(0.9, 0.8, 0.1),    # Yellow
            "SFR" => RGBf(0.1, 0.7, 0.7),    # Teal
            "MSR" => RGBf(0.6, 0.6, 0.6)     # Gray
        )

        for reactor_name in available_reactors
            pj_idx = findfirst(p -> p.name == reactor_name, pjs)
            if !isnothing(pj_idx)
                push!(reactor_types_in_scale, pjs[pj_idx].type)
            else
                push!(reactor_types_in_scale, "Unknown")
            end
        end

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

        # Create a separate axis for the reactor type indicator bars in the white space
        ax_type_indicator = Axis(gl[1, 3],
                                 yticks = (1:length(yticks_scale), fill("", length(yticks_scale))),
                                 limits = (0, 1, 0.5, length(yticks_scale) + 0.5))

        # Hide all decorations except the colored bars
        hidedecorations!(ax_type_indicator)
        hidespines!(ax_type_indicator)

        # Draw colored indicator bars for each reactor type
        for (j, rtype) in enumerate(reactor_types_in_scale)
            bar_color = get(type_indicator_colors_map, rtype, RGBf(0.8, 0.8, 0.8))
            # Draw a thick colored bar spanning the width of this axis
            poly!(ax_type_indicator,
                  Point2f[(0, j-0.4), (1, j-0.4), (1, j+0.4), (0, j+0.4)],
                  color = bar_color,
                  strokewidth = 0)
        end

        # Add type labels next to the indicator bars
        type_label_colors_map = Dict(
            "PWR" => RGBf(0.9, 0.5, 0.1),    # Orange
            "BWR" => RGBf(0.9, 0.5, 0.1),    # Orange (same as PWR)
            "HTR" => RGBf(0.9, 0.8, 0.1),    # Yellow
            "SFR" => RGBf(0.1, 0.7, 0.7),    # Teal
            "MSR" => RGBf(0.6, 0.6, 0.6)     # Gray
        )

        # Group consecutive reactors of the same type and add a single label
        type_groups = []
        if !isempty(reactor_types_in_scale)
            current_type = reactor_types_in_scale[1]
            start_idx = 1

            for j in 2:length(reactor_types_in_scale)
                if reactor_types_in_scale[j] != current_type
                    # Store the group info
                    mid_point = (start_idx + j - 1) / 2
                    push!(type_groups, (current_type, mid_point))
                    current_type = reactor_types_in_scale[j]
                    start_idx = j
                end
            end
            # Add the last group
            mid_point = (start_idx + length(reactor_types_in_scale)) / 2
            push!(type_groups, (current_type, mid_point))
        end

        # Draw type labels to the right of the indicator bars
        for (rtype, y_pos) in type_groups
            label_color = get(type_label_colors_map, rtype, :black)
            text!(ax_type_indicator, 1.2, y_pos,
                 text = rtype,
                 align = (:left, :center),
                 color = label_color,
                 fontsize = 10,
                 font = "Noto Sans Bold")
        end

        # Set column widths: heatmap gets most space, type indicator is narrow
        colsize!(gl, 3, Fixed(40))  # Type indicator column fixed at 40 pixels

        # Add colorbar inside the GridLayout
        Colorbar(gl[1, 2], hmap_shapley, width = 12)
        colsize!(gl, 2, Auto(12))  # Keep colorbar slim
    end

    # Overall title
    Label(fig[0, :], title, fontsize = 20, font = "Noto Sans Bold", color = (:black, 0.25))

    # Add legend explaining the colored type indicator bars at the bottom
    legend_layout = fig[2, :] = GridLayout()

    # Create legend elements for reactor types (using PolyElement for colored bars)
    legend_elements = []
    legend_labels = []

    # Define indicator colors for legend (matching the plot)
    legend_indicator_colors = [
        ("PWR/BWR", RGBf(0.9, 0.5, 0.1)),  # Orange
        ("HTR", RGBf(0.9, 0.8, 0.1)),      # Yellow
        ("SFR", RGBf(0.1, 0.7, 0.7))       # Teal
    ]

    for (type_name, indicator_color) in legend_indicator_colors
        push!(legend_elements, PolyElement(color = indicator_color, strokewidth = 0))
        push!(legend_labels, type_name)
    end

    Legend(legend_layout[1, 1],
           legend_elements,
           legend_labels,
           "Reactor Type Indicators (colored bars on the right side of each panel)",
           orientation = :horizontal,
           tellwidth = false,
           tellheight = true,
           framevisible = true,
           labelsize = 11,
           titlesize = 12)

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
                  xlabel = "[EUR2025/MW]",
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
              xlabel = "LCOE Threshold [EUR2025/MWh]",
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
                 ylabel="LCOE [EUR2025/MWh]",
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
                   label="Regional Mean: $(round(region_mean, digits=1)) EUR2025/MWh")

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
             ylabel="LCOE [EUR2025/MWh]",
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
- Figure with WACC (%) vs median LCOE (EUR2025/MWh) for Micro, SMR, and Large reactors
"""
function wacc_sensitivity_plot(outputpath, opt_scaling, pjs_dat, wacc_bin_centers; show_ci::Bool=false)

    @info("Generating WACC sensitivity plot from existing results (0 new simulations)")
    if show_ci
        @info("  Including 10th-90th percentile confidence intervals")
    end

    # Read saved results
    lcoe_results = CSV.read("$outputpath/mcs-lcoe_results-$opt_scaling.csv", DataFrame)
    wacc_values = CSV.read("$outputpath/mcs-wacc_values-$opt_scaling.csv", DataFrame)

    # Storage for binned results
    wacc_percentages = Float64[]
    micro_median_lcoe = Float64[]
    smr_median_lcoe = Float64[]
    large_median_lcoe = Float64[]

    # NEW: Storage for confidence intervals
    micro_p10_lcoe = Float64[]
    micro_p90_lcoe = Float64[]
    smr_p10_lcoe = Float64[]
    smr_p90_lcoe = Float64[]
    large_p10_lcoe = Float64[]
    large_p90_lcoe = Float64[]

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

        # NEW: Calculate confidence intervals
        push!(micro_p10_lcoe, isempty(micro_lcoe_bin) ? NaN : quantile(micro_lcoe_bin, 0.10))
        push!(micro_p90_lcoe, isempty(micro_lcoe_bin) ? NaN : quantile(micro_lcoe_bin, 0.90))
        push!(smr_p10_lcoe, isempty(smr_lcoe_bin) ? NaN : quantile(smr_lcoe_bin, 0.10))
        push!(smr_p90_lcoe, isempty(smr_lcoe_bin) ? NaN : quantile(smr_lcoe_bin, 0.90))
        push!(large_p10_lcoe, isempty(large_lcoe_bin) ? NaN : quantile(large_lcoe_bin, 0.10))
        push!(large_p90_lcoe, isempty(large_lcoe_bin) ? NaN : quantile(large_lcoe_bin, 0.90))

        @info("    Bin counts: Micro=$(length(micro_lcoe_bin)), SMR=$(length(smr_lcoe_bin)), Large=$(length(large_lcoe_bin))")
    end

    # Create plot with consistent styling
    fig = Figure(size=(1000, 700))
    ax = Axis(fig[1,1],
             xlabel="Discount Rate (WACC) [%]",
             ylabel="Median LCOE [EUR2025/MWh]",
             title="LCOE Sensitivity to Discount Rate by Reactor Scale ($opt_scaling scaling)")

    # Use consistent color scheme matching learning curves and other plots
    scale_colors = Dict(
        "Micro" => RGBf(0.2, 0.4, 0.8),    # Blue
        "SMR" => RGBf(0.8, 0.4, 0.0),      # Orange
        "Large" => RGBf(0.0, 0.6, 0.4)     # Teal/Green
    )

    # NEW: Plot confidence bands first (if requested)
    if show_ci
        band!(ax, wacc_percentages, large_p10_lcoe, large_p90_lcoe,
              color = (scale_colors["Large"], 0.2))
        band!(ax, wacc_percentages, smr_p10_lcoe, smr_p90_lcoe,
              color = (scale_colors["SMR"], 0.2))
        band!(ax, wacc_percentages, micro_p10_lcoe, micro_p90_lcoe,
              color = (scale_colors["Micro"], 0.2))
    end

    # Plot lines for each scale with thicker lines
    lines!(ax, wacc_percentages, micro_median_lcoe,
          label="Micro Reactors", linewidth=3.5, color=scale_colors["Micro"])
    lines!(ax, wacc_percentages, smr_median_lcoe,
          label="SMR Reactors", linewidth=3.5, color=scale_colors["SMR"])
    lines!(ax, wacc_percentages, large_median_lcoe,
          label="Large Reactors", linewidth=3.5, color=scale_colors["Large"])

    # Add markers with white stroke for visibility
    scatter!(ax, wacc_percentages, micro_median_lcoe,
            color=scale_colors["Micro"], markersize=12,
            strokewidth=1.5, strokecolor=:white)
    scatter!(ax, wacc_percentages, smr_median_lcoe,
            color=scale_colors["SMR"], markersize=12,
            strokewidth=1.5, strokecolor=:white)
    scatter!(ax, wacc_percentages, large_median_lcoe,
            color=scale_colors["Large"], markersize=12,
            strokewidth=1.5, strokecolor=:white)

    # Add legend with better positioning
    Legend(fig[1,2], ax, "Reactor Scale",
          framevisible=true, labelsize=12, titlesize=13)

    # Grid for readability
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridstyle = :dash
    ax.ygridstyle = :dash
    ax.xgridcolor = (:gray, 0.3)
    ax.ygridcolor = (:gray, 0.3)

    # Add more granular y-axis ticks for better readability
    # Determine appropriate range based on actual data
    all_values = vcat(micro_median_lcoe, smr_median_lcoe, large_median_lcoe)
    all_values = all_values[.!isnan.(all_values)]  # Remove NaNs

    y_min = floor(minimum(all_values) / 50) * 50  # Round down to nearest 50
    y_max = ceil(maximum(all_values) / 50) * 50   # Round up to nearest 50

    # Create tick marks every 50 EUR/MWh for granularity
    ytick_values = collect(y_min:50:y_max)
    ax.yticks = ytick_values

    @info("  Y-axis range: $(y_min) - $(y_max) EUR2025/MWh")
    @info("  Y-axis ticks: $(length(ytick_values)) marks at 50 EUR/MWh intervals")
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
- Y-axis: LCOE [EUR2025/MWh]
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
                ylabel = "LCOE [EUR2025/MWh]",
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
- `thresholds`: Vector of LCOE thresholds to evaluate (default: 0-300 EUR2025/MWh)

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
              xlabel = "LCOE Threshold [EUR2025/MWh]",
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

"""
    lcoe_comparison_from_mcs(lcoe_results::DataFrame, pjs::Vector; opt_scaling::String="rothwell")

Create LCOE comparison plot sourced directly from Monte Carlo simulation results.

This replaces the deterministic calculation approach with percentiles from the 
full MCS that includes:
- Investment cost uncertainty (including Large reactor overrun multiplier)
- WACC uncertainty [0.04, 0.10]
- Construction time uncertainty (scale-specific ranges)
- Load factor uncertainty (reactor-type-specific ranges)

# Arguments
- `lcoe_results`: DataFrame with LCOE MCS results (columns = reactor names)
- `pjs`: Vector of project objects containing reactor metadata
- `opt_scaling`: Scaling method name for labeling

# Returns
- Figure with percentile ranges for each reactor overlaid with Lazard data

Creates horizontal range plot showing:
- P10-P90 range (light gray, wide band)
- P25-P75 range (IQR, darker gray band)
- P50 median (black vertical line marker)
- Grouped by scale (Micro, SMR, Large)
- Overlaid with Lazard renewable/conventional energy LCOE ranges
"""
function lcoe_comparison_from_mcs(lcoe_results::DataFrame, pjs::Vector;
                                   opt_scaling::String="rothwell")
    
    # Load external LCOE data (Lazard LCOE 18.0, 2025 edition)
    inputpath = "_input"
    lcoe_dat = CSV.read("$inputpath/lcoe_data.csv", DataFrame)
    
    # Group reactors by scale
    scale_groups = Dict("Micro" => String[], "SMR" => String[], "Large" => String[])
    
    for pj in pjs
        if haskey(scale_groups, pj.scale) && (pj.name in names(lcoe_results))
            push!(scale_groups[pj.scale], pj.name)
        end
    end
    
    # Calculate percentiles for each reactor
    lcoe_stats = DataFrame(
        reactor = String[],
        scale = String[],
        p10 = Float64[],
        p25 = Float64[],
        p50 = Float64[],
        p75 = Float64[],
        p90 = Float64[]
    )
    
    scale_order = ["Micro", "SMR", "Large"]
    
    for scale in scale_order
        reactors = get(scale_groups, scale, String[])
        for reactor in reactors
            data = lcoe_results[!, reactor]
            push!(lcoe_stats, (
                reactor = reactor,
                scale = scale,
                p10 = quantile(data, 0.10),
                p25 = quantile(data, 0.25),
                p50 = quantile(data, 0.50),
                p75 = quantile(data, 0.75),
                p90 = quantile(data, 0.90)
            ))
        end
    end
    
    # Create figure
    fig = Figure(size=(1400, 1000))
    ax = Axis(fig[1,1],
             xlabel = "LCOE [EUR2025/MWh]",
             xscale = log10,
             ylabel = "Technology",
             xlabelsize = 16,
             ylabelsize = 16)
    
    xlims!(ax, 10, 25000)
    
    # Build y-axis positions and labels
    y_pos = 1
    ytick_labels = String[]
    ytick_positions = Int[]
    
    # Color scheme by scale (same as other plots)
    scale_colors = Dict(
        "Micro" => (1.0, 0.42, 0.42),  # Orange/coral
        "SMR" => (0.31, 0.80, 0.77),   # Cyan/blue
        "Large" => (0.58, 0.88, 0.83)  # Light cyan/teal
    )
    
    # Plot reactor data (grouped by scale)
    for scale in ["Large", "SMR", "Micro"]  # Reverse order (top to bottom)
        scale_data = filter(row -> row.scale == scale, lcoe_stats)
        
        if nrow(scale_data) == 0
            continue
        end
        
        scale_color = get(scale_colors, scale, (:gray, 0.5))
        
        for row in eachrow(scale_data)
            # P10-P90 range (wide transparent band)
            rangebars!(ax, [y_pos], [row.p10], [row.p90],
                      color=(scale_color, 0.3), linewidth=16,
                      whiskerwidth=10, direction=:x)
            
            # P25-P75 range (IQR, darker)
            rangebars!(ax, [y_pos], [row.p25], [row.p75],
                      color=(scale_color, 0.6), linewidth=16,
                      whiskerwidth=10, direction=:x)
            
            # Median (vertical line marker)
            scatter!(ax, [row.p50], [y_pos],
                    color=:black, marker=:vline,
                    markersize=24, strokewidth=2)
            
            push!(ytick_labels, row.reactor)
            push!(ytick_positions, y_pos)
            y_pos += 1
        end

        # Add spacing between scales (skip a row)
        y_pos += 1
    end

    # Add Lazard renewable/conventional energy ranges
    for row in eachrow(lcoe_dat)
        rangebars!(ax, [y_pos], [row.lower_bound], [row.upper_bound],
                  color=(:darkblue, 0.5), linewidth=16,
                  whiskerwidth=10, direction=:x)

        push!(ytick_labels, row.technology)
        push!(ytick_positions, y_pos)
        y_pos += 1
    end
    
    ax.yticks = (ytick_positions, ytick_labels)
    ax.yticklabelsize = 11
    ax.xgridvisible = true
    ax.ygridvisible = false
    
    # Add dividing lines
    n_reactors = nrow(lcoe_stats)
    hlines!(ax, [n_reactors + 0.75], color=:red, linestyle=:dash, linewidth=2)
    
    # Add title and caption
    Label(fig[0, :], "LCOE Comparison: Monte Carlo Simulation vs Benchmarks ($opt_scaling scaling)",
          fontsize=18, font="Noto Sans Bold")
    
    Label(fig[2, :], "Nuclear bands: 10th-90th (light) and 25th-75th (dark) percentiles from MCS (1M samples). Black line = median. Benchmark data: Lazard LCOE 18.0 (2025).",
          fontsize=11)
    
    return fig
end

"""
    lcoe_comparison_horizontal(lcoe_results::DataFrame, pjs::Vector, opt_scaling::String)

Create LCOE comparison plot matching the original fig_lcoe_comparison style.
Extends the original plot to include Large and Micro reactors above SMRs.

Uses rangebars to show Q25-Q75 ranges, exactly like the original code.
Includes Lazard 2025 benchmark data for renewables and conventionals.

# Arguments
- `lcoe_results::DataFrame`: Monte Carlo LCOE results (columns = reactors)
- `pjs::Vector`: Vector of project objects with metadata
- `opt_scaling::String`: Scaling method name for labeling

# Returns
- `Figure`: Makie figure object matching original style
"""
function lcoe_comparison_horizontal(lcoe_results::DataFrame, pjs::Vector, opt_scaling::String)
    # NOTE: lcoe_results parameter kept for backwards compatibility but no longer used
    # Now reads from pre-calculated summary statistics for consistency

    # Read Lazard LCOE data from CSV (same as original code)
    inputpath = "_input"
    outputpath = "_output"
    lcoe_dat = CSV.read("$inputpath/lcoe_data.csv", DataFrame)

    # FIXED: Read pre-calculated summary statistics instead of recalculating from raw data
    # This ensures consistency with published summary statistics
    lcoe_summary = CSV.read("$outputpath/mcs-lcoe_summary-$opt_scaling.csv", DataFrame)

    # Create a lookup dictionary for Q25/Q75 by reactor name
    lcoe_bounds_dict = Dict{String, Tuple{Float64, Float64}}()
    for row in eachrow(lcoe_summary)
        lcoe_bounds_dict[row.variable] = (row.q25, row.q75)
    end

    # FIXED: Calculate aggregated ranges BY TYPE for each scale using reactor names
    # Previous code had broken indexing - reactors aren't ordered by scale in the file

    # Helper function to get bounds for reactors matching scale and types
    function get_bounds_for_group(scale::String, types::Vector{String})
        matching_names = [pj.name for pj in pjs if pj.scale == scale && pj.type in types]
        if isempty(matching_names)
            return (0.0, 0.0)
        end
        q25_vals = [lcoe_bounds_dict[name][1] for name in matching_names if haskey(lcoe_bounds_dict, name)]
        q75_vals = [lcoe_bounds_dict[name][2] for name in matching_names if haskey(lcoe_bounds_dict, name)]
        return (minimum(q25_vals), maximum(q75_vals))
    end

    # LARGE reactors by type
    large_pwr_lower, large_pwr_upper = get_bounds_for_group("Large", ["PWR", "BWR"])
    large_htr_lower, large_htr_upper = get_bounds_for_group("Large", ["HTR"])
    large_sfr_lower, large_sfr_upper = get_bounds_for_group("Large", ["SFR"])

    # MICRO reactors by type
    micro_pwr_lower, micro_pwr_upper = get_bounds_for_group("Micro", ["PWR", "BWR"])
    micro_htr_lower, micro_htr_upper = get_bounds_for_group("Micro", ["HTR"])
    micro_sfr_lower, micro_sfr_upper = get_bounds_for_group("Micro", ["SFR"])

    # SMR reactors by type
    smr_pwr_lower, smr_pwr_upper = get_bounds_for_group("SMR", ["PWR", "BWR"])
    smr_htr_lower, smr_htr_upper = get_bounds_for_group("SMR", ["HTR"])
    smr_sfr_lower, smr_sfr_upper = get_bounds_for_group("SMR", ["SFR"])

    # Build plot data with reactor types for each scale
    lcoe_plot_data = vcat(
        select(lcoe_dat, [:technology, :lower_bound, :upper_bound]),
        DataFrame(
            technology = [
                "BWR & PWR Large", "HTR Large", "SFR Large",
                "BWR & PWR Micro", "HTR Micro", "SFR Micro",
                "BWR & PWR SMRs", "HTR SMRs", "SFR SMRs"
            ],
            lower_bound = [
                large_pwr_lower, large_htr_lower, large_sfr_lower,
                micro_pwr_lower, micro_htr_lower, micro_sfr_lower,
                smr_pwr_lower, smr_htr_lower, smr_sfr_lower
            ],
            upper_bound = [
                large_pwr_upper, large_htr_upper, large_sfr_upper,
                micro_pwr_upper, micro_htr_upper, micro_sfr_upper,
                smr_pwr_upper, smr_htr_upper, smr_sfr_upper
            ]
        )
    )

    # Define LCOE plot (EXACTLY like original code)
    if opt_scaling == "manufacturer"
        plot_scaling = "Manufacturer"
    elseif opt_scaling == "roulstone"
        plot_scaling = "Roulstone"
    elseif opt_scaling == "rothwell"
        plot_scaling = "Rothwell"
    elseif opt_scaling == "uniform"
        plot_scaling = "uniform"
    else
        @error("scaling not defined")
        plot_scaling = opt_scaling
    end

    xlabel = "[EUR2025/MWh]"
    yticks = lcoe_plot_data[!, :technology]

    # Color assignment - now with 3 types per scale
    n_renewables = 7  # From Lazard CSV
    n_conventionals = 4  # From Lazard CSV
    col = vcat(
        fill(1, n_renewables),      # Renewables (green)
        fill(2, n_conventionals),    # Conventionals (blue)
        3, 4, 5,  # Large: PWR/BWR, HTR, SFR
        6, 7, 8,  # Micro: PWR/BWR, HTR, SFR
        9, 10, 11  # SMR: PWR/BWR, HTR, SFR
    )

    colormap = [:darkgreen, :darkblue]

    fig_lcoe_comparison = Figure()
    ax_lcoe = Axis(fig_lcoe_comparison[1,1],
                   yticks = (1:length(yticks), yticks),
                   xscale = log10,
                   xlabel = xlabel)

    xlims!(10, 25000)

    # Rangebars (EXACTLY like original)
    rangebars!(ax_lcoe, 1:length(yticks), lcoe_plot_data[!, 2], lcoe_plot_data[!, 3],
              linewidth = 6, whiskerwidth = 8, direction = :x, color = col)

    # Separator lines - add lines between each scale group
    n_large = 3  # 3 types
    n_micro = 3  # 3 types
    n_smr = 3    # 3 types

    separators = [
        n_renewables + 0.5,  # After renewables (red)
        n_renewables + n_conventionals + 0.5,  # After conventionals (red)
        n_renewables + n_conventionals + n_large + 0.5,  # After Large (orange)
        n_renewables + n_conventionals + n_large + n_micro + 0.5  # After Micro (orange)
    ]

    line_colors = [:red, :red, :orange, :orange]

    hlines!(ax_lcoe, separators, linestyle = :dash, color = line_colors, linewidth = 1.5)

    # Text labels for categories
    renewables_center = n_renewables / 2
    conventionals_center = n_renewables + (n_conventionals / 2) + 0.5
    large_center = n_renewables + n_conventionals + (n_large / 2) + 0.5
    micro_center = n_renewables + n_conventionals + n_large + (n_micro / 2) + 0.5
    smr_center = n_renewables + n_conventionals + n_large + n_micro + (n_smr / 2) + 0.5

    text!([15000, 15000, 16, 16, 16],
          [renewables_center, conventionals_center, large_center, micro_center, smr_center];
          text = ["Renewables\n(LAZARD 2025)", "Conventionals\n(LAZARD 2025)",
                  "Large", "Micro", "SMR"],
          align = (:center, :center),
          justification = :center,
          rotation = π/2)

    # Value labels (EXACTLY like original)
    text!(lcoe_plot_data[!, 2], 1:length(yticks),
         text = string.(round.(Int, lcoe_plot_data[!, 2])),
         align = (:right, :center),
         offset = (-10, 0))

    text!(lcoe_plot_data[!, 3], 1:length(yticks),
         text = string.(round.(Int, lcoe_plot_data[!, 3])),
         align = (:left, :center),
         offset = (10, 0))

    # Title (EXACTLY like original)
    Label(fig_lcoe_comparison[1, 1, Top()], "LCOE Comparison",
         font = "Noto Sans Bold", padding = (0, 6, 6, 0))

    return fig_lcoe_comparison
end

"""
    shapley_heatmap_threepanel(shapley_results::DataFrame, pjs::Vector)

Create three-panel heatmap showing Shapley sensitivity indices by reactor scale.

# Features
- Three panels: Micro | SMR | Large
- Rows: wacc, construction_time, loadfactor, investment
- Color scale: 0.0 (white) to 1.0 (dark blue)
- Numerical values displayed in cells (2 decimals)
- Verifies efficiency property: sum ≈ 1.0

# Arguments
- `shapley_results::DataFrame`: Shapley indices with columns for each reactor
- `pjs::Vector`: Vector of project objects

# Returns
- `Figure`: Makie figure object
"""
function shapley_heatmap_threepanel(shapley_results::DataFrame, pjs::Vector)
    
    # Group reactors by scale
    scales = ["Micro", "SMR", "Large"]
    fig = Figure(size=(1400, 600))
    
    for (panel_idx, scale) in enumerate(scales)
        # Filter reactors by scale
        scale_reactors = [pj for pj in pjs if pj.scale == scale]
        
        if isempty(scale_reactors)
            continue
        end
        
        # Extract Shapley values for this scale
        n_reactors = length(scale_reactors)
        n_vars = 4  # wacc, construction_time, loadfactor, investment
        
        shapley_matrix = zeros(n_vars, n_reactors)
        reactor_names = String[]
        
        for (j, pj) in enumerate(scale_reactors)
            push!(reactor_names, pj.name)
            col_name = pj.name
            if col_name in names(shapley_results)
                shapley_matrix[1, j] = shapley_results[findfirst(==("wacc"), shapley_results.var), col_name]
                shapley_matrix[2, j] = shapley_results[findfirst(==("construction_time"), shapley_results.var), col_name]
                shapley_matrix[3, j] = shapley_results[findfirst(==("capacity factor"), shapley_results.var), col_name]
                shapley_matrix[4, j] = shapley_results[findfirst(==("scaling"), shapley_results.var), col_name]
            end
        end
        
        # Create heatmap
        ax = Axis(fig[1, panel_idx],
                 title = "$scale Reactors",
                 xlabel = "",
                 ylabel = panel_idx == 1 ? "Variable" : "",
                 xticks = (1:n_reactors, reactor_names),
                 xticklabelrotation = π/4,
                 yticks = (1:4, ["WACC", "Construction Time", "Load Factor", "Investment"]))
        
        # Use same colormap as Sobol sensitivity plots (matching thesis style)
        # Gradient from light (low values) to dark (high values)
        hm = heatmap!(ax, 1:n_reactors, 1:n_vars, shapley_matrix',
                     colormap = :RdYlBu_9,  # Red-Yellow-Blue gradient matching Sobol style
                     colorrange = (0.0, 1.0))

        # Add numerical values
        for i in 1:n_reactors
            for j in 1:n_vars
                val = shapley_matrix[j, i]
                text!(ax, i, j,
                     text = string(round(val, digits=2)),
                     align = (:center, :center),
                     fontsize = 10,
                     color = val > 0.4 ? :white : :black)  # Adjusted threshold for better contrast
            end
        end
        
        # Add colorbar to rightmost panel
        if panel_idx == length(scales)
            Colorbar(fig[1, panel_idx+1], hm, label = "Shapley Index")
        end
        
        # Verify efficiency property
        for i in 1:n_reactors
            col_sum = sum(shapley_matrix[:, i])
            if abs(col_sum - 1.0) > 0.05
                @warn "Shapley indices for $(reactor_names[i]) sum to $col_sum (expected ≈1.0)"
            end
        end
    end
    
    Label(fig[0, :], "Shapley Sensitivity Indices by Reactor Scale",
          fontsize = 18, font = "Noto Sans Bold")
    
    return fig
end

"""
    learning_curves_smr(pjs::Vector, wacc_range::Vector, electricity_price_mean::Float64, 
                        opt_scaling::String, outputpath::String)

Generate learning curves for SMR reactors showing LCOE reduction with cumulative units.

# Features
- Separate subplot for each SMR reactor
- Three curves: Pessimistic (5% LR), Base (10% LR), Optimistic (15% LR)
- X-axis: Cumulative Units (1-30)
- Y-axis: LCOE [EUR2025/MWh]
- Shaded confidence bands (P10-P90)
- Color-coded: Pessimistic=red, Base=orange, Optimistic=green

# Mathematical Basis
Wright's Law: C_NOAK = C_FOAK × (1 - LR)^(log₂(N))

# Arguments
- `pjs::Vector`: Vector of project objects
- `wacc_range::Vector`: [min, max] WACC range
- `electricity_price_mean::Float64`: Mean electricity price
- `opt_scaling::String`: Scaling method
- `outputpath::String`: Path to read baseline LCOE data

# Returns
- `Figure`: Makie figure object
"""
function learning_curves_smr(pjs::Vector, wacc_range::Vector, electricity_price_mean::Float64,
                              opt_scaling::String, outputpath::String)

    # Use provided reactor list (can be filtered by scale)
    all_reactors = pjs

    if isempty(all_reactors)
        @warn "No reactors found"
        return Figure()
    end

    # Determine scale for title
    scales_in_data = unique([pj.scale for pj in all_reactors])
    scale_label = length(scales_in_data) == 1 ? scales_in_data[1] : "All"

    # Learning parameters
    N_values = [1, 2, 4, 6, 8, 12, 16, 20, 24, 30]
    learning_rates = [0.05, 0.10, 0.15]  # Pessimistic, Base, Optimistic

    # Line widths: baseline (10%) thick, others thinner
    lr_linewidths = Dict(0.05 => 2.0, 0.10 => 3.5, 0.15 => 2.0)
    lr_labels = Dict(0.05 => "LR 5%", 0.10 => "LR 10% (Base)", 0.15 => "LR 15%")

    # Updated color scheme - darker colors as requested
    reactor_type_colors = Dict(
        "BWR" => RGBf(0.8, 0.4, 0.0),    # Darker orange
        "PWR" => RGBf(0.8, 0.4, 0.0),    # Darker orange (same as BWR)
        "HTR" => RGBf(0.7, 0.6, 0.0),    # Darker yellow
        "SFR" => RGBf(0.0, 0.5, 0.5),    # Darker teal
        "MSR" => RGBf(0.5, 0.5, 0.5)     # Gray
    )

    # Create grid layout - 2 columns × 5 rows (PORTRAIT orientation for PDF appendix)
    # This gives a tall figure that fits better on A4/Letter pages
    n_reactors = length(all_reactors)
    n_cols = 2   # 2 columns wide
    n_rows = Int(ceil(n_reactors / n_cols))  # 5 rows for 10 reactors

    # Figure size: 2 columns × 5 rows = portrait orientation
    # Width = 500 per column × 2 = 1000px
    # Height = 400 per row × 5 = 2000px (tall portrait)
    fig = Figure(size=(500 * n_cols, 400 * n_rows))

    # Try to load baseline LCOE data
    baseline_file = "$outputpath/mcs-lcoe_summary-$opt_scaling.csv"
    baseline_lcoe = Dict{String, Float64}()
    baseline_std = Dict{String, Float64}()

    if isfile(baseline_file)
        summary_data = CSV.read(baseline_file, DataFrame)
        for (i, pj) in enumerate(all_reactors)
            if "variable" in names(summary_data)
                idx = findfirst(==(pj.name), summary_data.variable)
                if !isnothing(idx)
                    if "median" in names(summary_data)
                        baseline_lcoe[pj.name] = summary_data[idx, "median"]
                    end
                    if "std" in names(summary_data)
                        baseline_std[pj.name] = summary_data[idx, "std"]
                    end
                end
            end
        end
    end

    # DISABLED: Vendor target LCOE feature (not used with new learning curve functions)
    # The new learning_curve_single_reactor() and learning_curves_grid_appendix()
    # calculate vendor baseline on-the-fly using calculate_vendor_baseline_lcoe_simple()
    vendor_target_lcoe = Dict{String, Float64}()

    # # Load vendor target LCOE (manufacturer-based baseline)
    # vendor_file = "$outputpath/mcs-lcoe_summary-manufacturer-baseline.csv"
    # if isfile(vendor_file)
    #     vendor_data = CSV.read(vendor_file, DataFrame)
    #     if "variable" in names(vendor_data) && "median" in names(vendor_data)
    #         for row in eachrow(vendor_data)
    #             vendor_target_lcoe[row.variable] = row.median
    #         end
    #     end
    # else
    #     @warn "Vendor target file not found: $vendor_file"
    # end

    for (idx, pj) in enumerate(all_reactors)
        row = div(idx - 1, n_cols) + 1
        col = mod(idx - 1, n_cols) + 1

        ax = Axis(fig[row, col],
                 title = pj.name,
                 xlabel = "Cumulative Units (N)",
                 ylabel = row == 1 && col == 1 ? "LCOE [EUR2025/MWh]" : "")

        # Get baseline FOAK LCOE and uncertainty
        foak_lcoe = get(baseline_lcoe, pj.name, 100.0)
        foak_std = get(baseline_std, pj.name, 20.0)  # Default std if not found

        # Get color for this reactor type
        reactor_color = get(reactor_type_colors, pj.type, :gray)

        # Calculate vendor baseline using OPTIMISTIC assumptions (7% WACC, 3yr construction)
        vendor_lcoe = calculate_vendor_baseline_lcoe_simple(pj, 0.07, 3)

        # Add vendor baseline as horizontal reference line
        hlines!(ax, [vendor_lcoe];
               color = :red,
               linestyle = :dashdot,
               linewidth = 2.0,
               label = idx == 1 ? "Vendor Baseline (optimistic)" : "")

        # Plot learning curves with confidence bands
        for LR in learning_rates
            lcoe_curve = Float64[]
            lcoe_lower = Float64[]  # 95% CI lower bound
            lcoe_upper = Float64[]  # 95% CI upper bound

            for N in N_values
                # Apply Wright's Law
                multiplier = (1 - LR)^(log2(N))
                noak_lcoe = foak_lcoe * multiplier
                push!(lcoe_curve, noak_lcoe)

                # Propagate uncertainty (simplified - assume std scales with mean)
                noak_std = foak_std * multiplier
                push!(lcoe_lower, noak_lcoe - 1.96 * noak_std)  # 95% CI
                push!(lcoe_upper, noak_lcoe + 1.96 * noak_std)
            end

            # Plot confidence band for ALL learning rates (as requested)
            band_alpha = LR == 0.10 ? 0.25 : 0.15  # Base case slightly more opaque
            band!(ax, N_values, lcoe_lower, lcoe_upper,
                 color = (reactor_color, band_alpha))

            # Plot line - all solid, varying thickness
            line_width = lr_linewidths[LR]
            lines!(ax, N_values, lcoe_curve,
                  color = reactor_color,
                  linestyle = :solid,  # All solid as requested
                  linewidth = line_width,
                  label = lr_labels[LR])

            # Add dots for each observation point
            scatter!(ax, N_values, lcoe_curve,
                    color = reactor_color,
                    markersize = LR == 0.10 ? 8 : 6,  # Larger dots for baseline
                    strokewidth = 1,
                    strokecolor = :white)

            # Mark learning threshold where curve meets vendor target
            if !isnothing(vendor_lcoe)
                idx_hit = findfirst(x -> x <= vendor_lcoe, lcoe_curve)
                if !isnothing(idx_hit)
                    N_hit = N_values[idx_hit]
                    vlines!(ax, [N_hit];
                           color = reactor_color,
                           linestyle = :dash,
                           linewidth = 1)
                    scatter!(ax, [N_hit], [lcoe_curve[idx_hit]];
                            color = reactor_color,
                            markersize = 7,
                            strokewidth = 1)
                end
            end
        end

        # Add legend to first subplot with reactor type color indicator
        if idx == 1
            axislegend(ax, position = :rt, framevisible = true,
                      labelsize = 10, titlesize = 11)
        end
    end

    # Add overall title with scale label
    Label(fig[0, :], "$scale_label Reactors: LCOE Reduction with Experience ($opt_scaling scaling)",
          fontsize = 18, font = "Noto Sans Bold")

    # Add global legend for reactor types
    legend_elements = []
    legend_labels = []
    for (rtype, color) in sort(collect(reactor_type_colors), by=x->x[1])
        # Check if this reactor type exists in our data
        if any(pj.type == rtype for pj in all_reactors)
            push!(legend_elements, LineElement(color = color, linewidth = 3))
            push!(legend_labels, rtype)
        end
    end

    if !isempty(legend_elements)
        Legend(fig[end+1, :], legend_elements, legend_labels,
              "Reactor Type Colors", orientation = :horizontal,
              framevisible = true, tellwidth = false, tellheight = true)
    end

    return fig
end

"""
    learning_curve_standalone(pj::project, wacc_range::Vector, electricity_price_mean::Float64,
                              opt_scaling::String, outputpath::String)

Create a standalone learning curve plot for a single reactor (for thesis main text).
Produces a larger, publication-quality figure with full legend and annotations.

# Arguments
- `pj::project`: Single reactor project object
- `wacc_range::Vector`: WACC range [min, max]
- `electricity_price_mean::Float64`: Mean electricity price
- `opt_scaling::String`: Scaling method identifier
- `outputpath::String`: Path to read baseline LCOE data

# Returns
- `Figure`: Makie figure object
"""
function learning_curve_standalone(pj::project, wacc_range::Vector, electricity_price_mean::Float64,
                                   opt_scaling::String, outputpath::String)

    @info "Creating standalone learning curve for $(pj.name)"

    # Learning parameters
    N_values = [1, 2, 4, 6, 8, 12, 16, 20, 24, 30]
    learning_rates = [0.05, 0.10, 0.15]  # Pessimistic, Base, Optimistic

    lr_linewidths = Dict(0.05 => 2.0, 0.10 => 3.5, 0.15 => 2.0)
    lr_labels = Dict(0.05 => "LR 5% (Pessimistic)", 0.10 => "LR 10% (Base)", 0.15 => "LR 15% (Optimistic)")

    # Color for reactor type
    reactor_type_colors = Dict(
        "BWR" => RGBf(0.8, 0.4, 0.0),
        "PWR" => RGBf(0.8, 0.4, 0.0),
        "HTR" => RGBf(0.7, 0.6, 0.0),
        "SFR" => RGBf(0.0, 0.5, 0.5),
        "MSR" => RGBf(0.5, 0.5, 0.5)
    )
    reactor_color = get(reactor_type_colors, pj.type, :gray)

    # Create single figure - 10% narrower, 5% less height for better PDF fit
    # Layout: title at top, plot in middle, legend at bottom
    fig = Figure(size=(920, 520))
    colsize!(fig.layout, 1, Relative(1))  # Force column to take full figure width

    # Load baseline LCOE from simulation results
    baseline_file = "$outputpath/mcs-lcoe_summary-$opt_scaling.csv"
    foak_lcoe = 100.0  # Default
    foak_std = 20.0    # Default

    if isfile(baseline_file)
        summary_data = CSV.read(baseline_file, DataFrame)
        if "variable" in names(summary_data)
            idx = findfirst(==(pj.name), summary_data.variable)
            if !isnothing(idx)
                if "median" in names(summary_data)
                    foak_lcoe = summary_data[idx, "median"]
                end
                if "std" in names(summary_data)
                    foak_std = summary_data[idx, "std"]
                end
            end
        end
    end

    # Calculate vendor baseline
    vendor_lcoe = calculate_vendor_baseline_lcoe_simple(pj, 0.07, 3)

    # Create axis with full labels - in row 2 to leave room for title
    # aspect = nothing disables any global theme aspect constraints
    # tellwidth/tellheight = true ensures axis participates in layout sizing
    ax = Axis(fig[2, 1],
             title = "$(pj.name): Learning Curve Analysis",
             titlesize = 16,
             xlabel = "Cumulative Units Built (N)",
             ylabel = "LCOE [EUR₂₀₂₅/MWh]",
             xlabelsize = 14,
             ylabelsize = 14,
             aspect = nothing,
             autolimitaspect = nothing,
             tellwidth = true,
             tellheight = true)

    # Store legend entries manually
    legend_entries = []
    legend_labels_list = []

    # Add vendor baseline as horizontal reference line
    hline = hlines!(ax, [vendor_lcoe];
           color = :red,
           linestyle = :dashdot,
           linewidth = 2.5)
    push!(legend_entries, [hline])
    push!(legend_labels_list, "Vendor Baseline (7% WACC, 90% CF)")

    # Plot learning curves with confidence bands
    for LR in learning_rates
        lcoe_curve = Float64[]
        lcoe_lower = Float64[]
        lcoe_upper = Float64[]

        for N in N_values
            # Apply Wright's Law
            multiplier = (1 - LR)^(log2(N))
            noak_lcoe = foak_lcoe * multiplier
            push!(lcoe_curve, noak_lcoe)

            # Propagate uncertainty
            noak_std = foak_std * multiplier
            push!(lcoe_lower, noak_lcoe - 1.96 * noak_std)
            push!(lcoe_upper, noak_lcoe + 1.96 * noak_std)
        end

        # Plot confidence band
        band_alpha = LR == 0.10 ? 0.25 : 0.15
        band!(ax, N_values, lcoe_lower, lcoe_upper,
             color = (reactor_color, band_alpha))

        # Plot line
        line = lines!(ax, N_values, lcoe_curve,
              color = reactor_color,
              linestyle = :solid,
              linewidth = lr_linewidths[LR])

        push!(legend_entries, [line])
        push!(legend_labels_list, lr_labels[LR])

        # Add dots for each point
        scatter!(ax, N_values, lcoe_curve,
                color = reactor_color,
                markersize = LR == 0.10 ? 10 : 7,
                strokewidth = 1,
                strokecolor = :white)

        # Mark crossing point with vendor baseline
        idx_hit = findfirst(x -> x <= vendor_lcoe, lcoe_curve)
        if !isnothing(idx_hit)
            N_hit = N_values[idx_hit]
            scatter!(ax, [N_hit], [lcoe_curve[idx_hit]];
                    color = :red,
                    markersize = 12,
                    marker = :star5,
                    strokewidth = 1,
                    strokecolor = :black)

            # Add annotation for crossing point (only for base case)
            if LR == 0.10
                text!(ax, N_hit + 1, lcoe_curve[idx_hit] + 10,
                      text = "Competitive at\nN = $N_hit units",
                      fontsize = 10,
                      color = :red,
                      align = (:left, :bottom))
            end
        end
    end

    # Add legend BELOW the plot (row 3) - larger text for PDF readability
    Legend(fig[3, 1], legend_entries, legend_labels_list,
           orientation = :horizontal,
           framevisible = true,
           labelsize = 14,
           patchsize = (20, 20),
           tellwidth = false,
           tellheight = true)

    # Add reactor type info as subtitle at TOP (row 1)
    # tellwidth = false prevents label from shrinking column width
    Label(fig[1, 1], "$(pj.type) | $(pj.scale) | $(Int(pj.plant_capacity)) MW",
          fontsize = 12, color = :gray50, tellwidth = false)

    # Adjust row sizes: small title, large plot, small legend
    rowsize!(fig.layout, 1, Auto(0.05))  # Title row
    rowsize!(fig.layout, 2, Auto(1.0))   # Plot row (main content)
    rowsize!(fig.layout, 3, Auto(0.1))   # Legend row

    return fig
end

"""
    threshold_probability_curves(lcoe_results::DataFrame, pjs::Vector)

Create S-curves showing probability of LCOE exceeding various thresholds.

# Features
- Smooth S-curves: P(LCOE > threshold) vs threshold
- X-axis: LCOE Threshold [EUR2025/MWh] from 0 to 300
- Y-axis: Probability of Exceeding [%] from 0 to 100
- Horizontal line at 50% probability
- Color-coded by reactor type
- All reactors on one plot

# Arguments
- `lcoe_results::DataFrame`: Monte Carlo LCOE results
- `pjs::Vector`: Vector of project objects

# Returns
- `Figure`: Makie figure object
"""
function threshold_probability_curves(lcoe_results::DataFrame, pjs::Vector)

    # Group reactors by scale
    scale_groups = Dict("Micro" => [], "SMR" => [], "Large" => [])
    for pj in pjs
        if haskey(scale_groups, pj.scale)
            push!(scale_groups[pj.scale], pj)
        end
    end

    # Color mapping by reactor type
    type_colors = Dict(
        "PWR" => RGBf(0.8, 0.2, 0.2),    # Red
        "BWR" => RGBf(0.9, 0.5, 0.0),    # Orange
        "HTR" => RGBf(0.9, 0.8, 0.0),    # Yellow
        "SFR" => RGBf(0.0, 0.6, 0.6),    # Cyan
        "MSR" => RGBf(0.6, 0.0, 0.6)     # Purple
    )

    # Threshold range (extended to 500 for better visualization)
    thresholds = 0.0:5.0:500.0

    # Create one figure per scale
    figures = Dict{String, Figure}()

    for scale in ["Micro", "SMR", "Large"]
        scale_pjs = scale_groups[scale]

        if isempty(scale_pjs)
            continue
        end

        fig = Figure(size=(1100, 750))
        ax = Axis(fig[1, 1],
                 xlabel = "LCOE Threshold [EUR₂₀₂₅/MWh]",
                 ylabel = "Probability of Exceeding Threshold [%]",
                 title = "$scale Reactors: LCOE Threshold Exceedance Probability",
                 titlesize = 22,
                 xlabelsize = 20,
                 ylabelsize = 20,
                 xticklabelsize = 18,
                 yticklabelsize = 18)

        # Plot each reactor in this scale
        for pj in scale_pjs
            lcoe_vals = lcoe_results[!, pj.name]
            prob_exceed = Float64[]

            for threshold in thresholds
                # Calculate empirical probability of exceeding threshold
                p_exceed = sum(lcoe_vals .> threshold) / length(lcoe_vals) * 100
                push!(prob_exceed, p_exceed)
            end

            # Get color for reactor type
            line_color = get(type_colors, pj.type, RGBf(0.5, 0.5, 0.5))

            lines!(ax, collect(thresholds), prob_exceed,
                  color = line_color,
                  linewidth = 3,
                  label = pj.name)
        end

        # Add 50% probability line
        hlines!(ax, [50.0], linestyle = :dash, color = :black, linewidth = 2)

        # Add legend with larger text
        Legend(fig[1, 2], ax, "Reactors",
               framevisible = true,
               labelsize = 16,
               titlesize = 18,
               patchsize = (30, 30))

        figures[scale] = fig
    end

    return figures
end

"""
    plot_lcoe_breakdown_by_scale(pjs::Vector, all_breakdown_results::Dict, opt_scaling::String)

Create stacked bar charts showing LCOE component breakdown by reactor scale.
Generates 3 figures (Micro, SMR, Large), each with stacked bars for reactor types within that scale.

# Arguments
- `pjs::Vector`: Vector of project objects
- `all_breakdown_results::Dict`: Dict mapping reactor names to breakdown results (from npv_lcoe with decompose=true)
- `opt_scaling::String`: Scaling method name for titles

# Returns
- Dict with keys "Micro", "SMR", "Large" containing Figure objects
"""
function plot_lcoe_breakdown_by_scale(pjs::Vector, all_breakdown_results::Dict, opt_scaling::String)

    @info "Creating LCOE breakdown stacked bar charts by scale..."

    # Group reactors by scale
    scale_groups = Dict("Micro" => [], "SMR" => [], "Large" => [])

    for pj in pjs
        if haskey(scale_groups, pj.scale) && haskey(all_breakdown_results, pj.name)
            push!(scale_groups[pj.scale], pj)
        end
    end

    # Color scheme for LCOE components (4-tier breakdown)
    # Consistent with thesis color scheme (matching learning curves and histograms)
    component_colors = Dict(
        "OCC" => RGBf(0.0, 0.5, 0.5),               # Teal (base capital cost)
        "IDC" => RGBf(0.8, 0.4, 0.0),               # Darker orange (financing cost - stands out)
        "Fixed O&M" => RGBf(0.7, 0.6, 0.0),         # Darker yellow (fixed operations)
        "Variable O&M + Fuel" => RGBf(0.8, 0.2, 0.2)  # Orangered (variable costs)
    )

    figures = Dict{String, Figure}()

    for (scale, reactors_in_scale) in scale_groups
        if isempty(reactors_in_scale)
            @info "  Skipping $scale (no reactors)"
            continue
        end

        @info "  Creating stacked bar chart for $scale (n=$(length(reactors_in_scale)) reactors)"

        # Collect data for this scale
        reactor_labels = String[]
        occ_values = Float64[]
        idc_values = Float64[]
        fixed_om_values = Float64[]
        variable_om_fuel_values = Float64[]

        for pj in reactors_in_scale
            res = all_breakdown_results[pj.name]

            # Use median values for the bars
            push!(reactor_labels, "$(pj.type)\n$(pj.name)")
            push!(occ_values, median(res.lcoe_occ))
            push!(idc_values, median(res.lcoe_idc))
            push!(fixed_om_values, median(res.lcoe_fixed_om))
            push!(variable_om_fuel_values, median(res.lcoe_variable_om_fuel))
        end

        # Create figure with 2-row grid layout for better readability
        n_reactors = length(reactor_labels)
        n_cols = min(5, n_reactors)  # Max 5 per row
        n_rows = ceil(Int, n_reactors / n_cols)

        # Figure width to match legend width - adjust this value to align plot with legend
        fig = Figure(size=(max(1100, 220 * n_cols), 500 * n_rows + 180), backgroundcolor=:white)

        # Force column to expand to full figure width
        colsize!(fig.layout, 1, Relative(1))

        # Title at top
        Label(fig[0, 1], "$scale Reactors: LCOE Component Breakdown ($opt_scaling scaling)",
              fontsize = 24, font = :bold, tellwidth = false)

        # Calculate max total for consistent y-axis across rows
        max_total = maximum(occ_values .+ idc_values .+ fixed_om_values .+ variable_om_fuel_values)

        # Create axes for each row
        for row in 1:n_rows
            start_idx = (row - 1) * n_cols + 1
            end_idx = min(row * n_cols, n_reactors)
            row_indices = start_idx:end_idx
            n_in_row = length(row_indices)

            ax = Axis(fig[row, 1],
                     ylabel = row == 1 ? "LCOE [EUR₂₀₂₅/MWh]" : "",
                     ylabelsize = 20,
                     yticklabelsize = 16,
                     xticks = (1:n_in_row, reactor_labels[row_indices]),
                     xticklabelrotation = π/4,
                     xticklabelsize = 16,
                     xgridvisible = false,
                     ygridvisible = true,
                     ygridcolor = (:gray, 0.2),
                     ygridstyle = :dash,
                     aspect = nothing,
                     autolimitaspect = nothing,
                     tellwidth = true,
                     tellheight = true)

            # Extract data for this row
            row_occ = occ_values[row_indices]
            row_idc = idc_values[row_indices]
            row_fixed_om = fixed_om_values[row_indices]
            row_variable = variable_om_fuel_values[row_indices]

            # Create stacked bars for this row
            barplot!(ax, 1:n_in_row, row_occ, color = component_colors["OCC"], width = 0.7)
            barplot!(ax, 1:n_in_row, row_idc, offset = row_occ, color = component_colors["IDC"], width = 0.7)
            barplot!(ax, 1:n_in_row, row_fixed_om, offset = row_occ .+ row_idc, color = component_colors["Fixed O&M"], width = 0.7)
            barplot!(ax, 1:n_in_row, row_variable, offset = row_occ .+ row_idc .+ row_fixed_om, color = component_colors["Variable O&M + Fuel"], width = 0.7)

            # Add simplified labels - only percentage inside bars (no EUR values to reduce clutter)
            for (local_i, global_i) in enumerate(row_indices)
                occ = occ_values[global_i]
                idc = idc_values[global_i]
                fixed_om = fixed_om_values[global_i]
                variable = variable_om_fuel_values[global_i]
                total = occ + idc + fixed_om + variable

                occ_pct = (occ / total) * 100
                idc_pct = (idc / total) * 100
                fixed_om_pct = (fixed_om / total) * 100
                variable_pct = (variable / total) * 100

                # Only show percentage labels if segment is large enough (>10%)
                if occ_pct > 10
                    text!(ax, local_i, occ / 2,
                         text = "$(round(Int, occ_pct))%",
                         align = (:center, :center), color = :white, fontsize = 14, font = :bold)
                end
                if idc_pct > 10
                    text!(ax, local_i, occ + idc / 2,
                         text = "$(round(Int, idc_pct))%",
                         align = (:center, :center), color = :black, fontsize = 14, font = :bold)
                end
                if fixed_om_pct > 10
                    text!(ax, local_i, occ + idc + fixed_om / 2,
                         text = "$(round(Int, fixed_om_pct))%",
                         align = (:center, :center), color = :black, fontsize = 14, font = :bold)
                end
                if variable_pct > 10
                    text!(ax, local_i, occ + idc + fixed_om + variable / 2,
                         text = "$(round(Int, variable_pct))%",
                         align = (:center, :center), color = :white, fontsize = 14, font = :bold)
                end

                # Total LCOE at top with full breakdown
                text!(ax, local_i, total + (total * 0.02),
                     text = "$(round(total, digits=0))",
                     align = (:center, :bottom),
                     color = :black,
                     fontsize = 16,
                     font = :bold)
            end

            # Set consistent y-axis limit for this row
            ylims!(ax, 0, max_total * 1.15)
        end

        # Add legend at the bottom (after all rows)
        legend_elements = [
            PolyElement(color = component_colors["OCC"], strokewidth = 0),
            PolyElement(color = component_colors["IDC"], strokewidth = 0),
            PolyElement(color = component_colors["Fixed O&M"], strokewidth = 0),
            PolyElement(color = component_colors["Variable O&M + Fuel"], strokewidth = 0)
        ]
        legend_labels = ["OCC (Overnight Construction Cost)", "IDC (Interest During Construction)", "Fixed O&M", "Variable O&M + Fuel"]

        Legend(fig[n_rows+1, 1],
               legend_elements,
               legend_labels,
               "LCOE Components",
               orientation = :horizontal,
               tellwidth = false,
               tellheight = true,
               framevisible = true,
               labelsize = 18,
               titlesize = 20,
               patchsize = (25, 25),
               padding = (15, 15, 10, 10))

        figures[scale] = fig
    end

    return figures
end

"""
    plot_idc_sensitivity(pj::project, wacc::Float64, base_construction_time::Int,
                         delayed_construction_time::Int)

Create waterfall chart comparing on-time vs delayed construction scenarios,
showing OCC and IDC components separately.

# Arguments
- `pj::project`: Reactor to analyze (e.g., BWRX-300)
- `wacc::Float64`: Discount rate (e.g., 0.07)
- `base_construction_time::Int`: On-time scenario (e.g., 5 years)
- `delayed_construction_time::Int`: Delayed scenario (e.g., 8 years)

# Returns
- Figure showing side-by-side comparison of capital cost breakdown
"""
function plot_idc_sensitivity(pj::project, wacc::Float64,
                               base_construction_time::Int,
                               delayed_construction_time::Int)

    # Use manufacturer OCC (nominal investment)
    occ = pj.investment * pj.plant_capacity  # Total overnight cost [EUR]

    # Calculate IDC for each scenario using explicit formula
    # IDC = OCC × [(r/2)×T_con + (r²/6)×T_con²]
    function calculate_idc(t_con, r, occ)
        idc_factor = (r/2) * t_con + (r^2/6) * t_con^2
        return occ * idc_factor
    end

    # Scenario 1: On-time construction
    idc_base = calculate_idc(base_construction_time, wacc, occ)
    total_base = occ + idc_base

    # Scenario 2: Delayed construction
    idc_delayed = calculate_idc(delayed_construction_time, wacc, occ)
    total_delayed = occ + idc_delayed

    # Calculate LCOE components
    # Annual generation
    loadfactor = mean(pj.loadfactor)
    annual_generation = pj.plant_capacity * loadfactor * 8760  # MWh/year
    lifetime = pj.time[2]  # Operating years

    # Capital recovery factor
    crf = wacc * (1 + wacc)^lifetime / ((1 + wacc)^lifetime - 1)

    # LCOE components [EUR/MWh]
    lcoe_occ_base = (occ * crf) / annual_generation
    lcoe_idc_base = (idc_base * crf) / annual_generation
    lcoe_capital_base = lcoe_occ_base + lcoe_idc_base

    lcoe_occ_delayed = (occ * crf) / annual_generation
    lcoe_idc_delayed = (idc_delayed * crf) / annual_generation
    lcoe_capital_delayed = lcoe_occ_delayed + lcoe_idc_delayed

    # Create figure with two side-by-side waterfalls
    fig = Figure(size=(1400, 700))

    # Scenario 1: On-time
    ax1 = Axis(fig[1,1],
              title = "On-Time Construction: $(base_construction_time) years",
              ylabel = "Capital LCOE [EUR2025/MWh]",
              xticks = (1:3, ["OCC", "IDC", "Total Capital"]),
              xticklabelrotation = 0)

    # Bars for on-time scenario
    barplot!(ax1, [1], [lcoe_occ_base], color=:steelblue3, width=0.6)
    barplot!(ax1, [2], [lcoe_idc_base], offset=[lcoe_occ_base],
             color=:coral2, width=0.6)
    barplot!(ax1, [3], [lcoe_capital_base], color=:darkblue, width=0.6)

    # Labels on-time
    text!(ax1, 1, lcoe_occ_base/2,
         text="$(round(lcoe_occ_base, digits=1))\n($(round(lcoe_occ_base/lcoe_capital_base*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=13, font=:bold)
    text!(ax1, 2, lcoe_occ_base + lcoe_idc_base/2,
         text="$(round(lcoe_idc_base, digits=1))\n($(round(lcoe_idc_base/lcoe_capital_base*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=13, font=:bold)
    text!(ax1, 3, lcoe_capital_base/2,
         text="$(round(lcoe_capital_base, digits=1))",
         align=(:center, :center), color=:white, fontsize=14, font=:bold)

    # Connecting line
    lines!(ax1, [1.3, 1.7], [lcoe_occ_base, lcoe_occ_base],
          color=:gray50, linestyle=:dash, linewidth=2)

    ylims!(ax1, 0, max(lcoe_capital_base, lcoe_capital_delayed) * 1.2)

    # Scenario 2: Delayed
    ax2 = Axis(fig[1,2],
              title = "Delayed Construction: $(delayed_construction_time) years (+$(delayed_construction_time - base_construction_time) years)",
              ylabel = "Capital LCOE [EUR2025/MWh]",
              xticks = (1:3, ["OCC", "IDC", "Total Capital"]),
              xticklabelrotation = 0)

    # Bars for delayed scenario
    barplot!(ax2, [1], [lcoe_occ_delayed], color=:steelblue3, width=0.6)
    barplot!(ax2, [2], [lcoe_idc_delayed], offset=[lcoe_occ_delayed],
             color=:orangered2, width=0.6)  # Darker red for emphasis
    barplot!(ax2, [3], [lcoe_capital_delayed], color=:darkred, width=0.6)

    # Labels delayed
    text!(ax2, 1, lcoe_occ_delayed/2,
         text="$(round(lcoe_occ_delayed, digits=1))\n($(round(lcoe_occ_delayed/lcoe_capital_delayed*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=13, font=:bold)
    text!(ax2, 2, lcoe_occ_delayed + lcoe_idc_delayed/2,
         text="$(round(lcoe_idc_delayed, digits=1))\n($(round(lcoe_idc_delayed/lcoe_capital_delayed*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=13, font=:bold)
    text!(ax2, 3, lcoe_capital_delayed/2,
         text="$(round(lcoe_capital_delayed, digits=1))",
         align=(:center, :center), color=:white, fontsize=14, font=:bold)

    # Connecting line
    lines!(ax2, [1.3, 1.7], [lcoe_occ_delayed, lcoe_occ_delayed],
          color=:gray50, linestyle=:dash, linewidth=2)

    ylims!(ax2, 0, max(lcoe_capital_base, lcoe_capital_delayed) * 1.2)

    # Add title and impact summary
    Label(fig[0, :],
          "Impact of Construction Delays on Capital Costs: $(pj.name) @ $(round(wacc*100, digits=0))% WACC",
          fontsize=20, font=:bold)

    # Calculate impact
    idc_increase = lcoe_idc_delayed - lcoe_idc_base
    total_increase = lcoe_capital_delayed - lcoe_capital_base
    pct_increase = (total_increase / lcoe_capital_base) * 100

    Label(fig[2, :],
          "Impact: +$(delayed_construction_time - base_construction_time) years delay → IDC increases by $(round(idc_increase, digits=1)) EUR/MWh → Total capital LCOE +$(round(total_increase, digits=1)) EUR/MWh (+$(round(pct_increase, digits=1))%)",
          fontsize=14, color=:darkred, font=:bold)

    # Legend
    legend_elements = [
        PolyElement(color=:steelblue3, strokewidth=0),
        PolyElement(color=:coral2, strokewidth=0),
        PolyElement(color=:orangered2, strokewidth=0)
    ]
    legend_labels = [
        "OCC (Overnight Construction Cost)",
        "IDC - On-time",
        "IDC - Delayed"
    ]

    Legend(fig[3, :],
           legend_elements,
           legend_labels,
           orientation=:horizontal,
           tellwidth=false,
           tellheight=true,
           framevisible=true)

    return fig
end

"""
    plot_idc_aggregate(pjs::Vector{project}, scale::String, wacc::Float64,
                       base_construction_time::Int, delayed_construction_time::Int)

Create single waterfall chart showing AGGREGATED capital cost breakdown for a reactor scale.

Shows average OCC for SMR/Large/Micro reactors and IDC under on-time vs delayed scenarios.

# Arguments
- `pjs::Vector{project}`: All reactor projects
- `scale::String`: "SMR", "Large", or "Micro"
- `wacc::Float64`: Discount rate (e.g., 0.07)
- `base_construction_time::Int`: On-time scenario (e.g., 5 years)
- `delayed_construction_time::Int`: Delayed scenario (e.g., 8 years)

# Returns
- Figure showing aggregated capital cost breakdown comparison
"""
function plot_idc_aggregate(pjs::Vector, scale::String, wacc::Float64,
                            base_construction_time::Int, delayed_construction_time::Int,
                            opt_scaling::String;
                            wacc_range::Vector=[0.03, 0.10],
                            electricity_price_mean::Float64=50.0,
                            construction_time_range::Vector{Int}=[3, 10])

    # Use BWRX-300 only to match standalone donut and aggregate breakdown
    scale_reactors = filter(p -> p.name == "BWRX-300", pjs)

    if isempty(scale_reactors)
        @error "BWRX-300 reactor not found"
        return nothing
    end

    p = scale_reactors[1]  # BWRX-300

    @info "Creating IDC sensitivity for BWRX-300 (using $opt_scaling scaling)"
    @info "  Using same simulation pipeline as donut/breakdown for consistent OCC"

    # ===== CONSISTENT WITH DONUT/BREAKDOWN PIPELINE =====
    # Run simulation to get capital costs (same as donut plot)
    n = 100_000  # Smaller sample for speed, still accurate for median
    all_capital = Float64[]
    all_wacc = Float64[]
    all_construction_time = Int[]

    rand_vars = gen_rand_vars(opt_scaling, n, wacc_range, electricity_price_mean, p;
                              construction_time_range=construction_time_range)
    disc_res = mc_run(n, p, rand_vars)
    res = npv_lcoe(disc_res, decompose=true)

    # Collect capital cost and parameters for OCC/IDC decomposition
    append!(all_capital, res.lcoe_capital)
    append!(all_wacc, rand_vars.wacc)
    append!(all_construction_time, rand_vars.construction_time)

    # Decompose capital into OCC and IDC using Wealer formula (same as donut)
    occ_illustrative, idc_illustrative = decompose_capital_to_occ_idc(
        all_capital, all_wacc, all_construction_time
    )

    # OCC0 = median illustrative OCC (this is our "OCC held constant" baseline)
    # This matches the donut/breakdown OCC value
    lcoe_occ = median(occ_illustrative)

    @info "  Median OCC (from simulation): $(round(lcoe_occ, digits=1)) EUR/MWh"
    @info "  Median IDC (from simulation): $(round(median(idc_illustrative), digits=1)) EUR/MWh"

    # ===== IDC SENSITIVITY CALCULATION =====
    # IDC factor formula (same as decompose_capital_to_occ_idc)
    idc_factor(T, r) = (r/2) * T + (r^2/6) * T^2

    # For sensitivity analysis, we need OCC in EUR (not EUR/MWh)
    # Back-calculate from LCOE using CRF and annual generation
    avg_capacity = p.plant_capacity
    avg_loadfactor = mean(p.loadfactor)
    avg_lifetime = p.time[2]
    annual_generation = avg_capacity * avg_loadfactor * 8760  # MWh/year
    crf = wacc * (1 + wacc)^avg_lifetime / ((1 + wacc)^avg_lifetime - 1)

    # OCC in EUR = lcoe_occ × annual_generation / crf
    occ_eur = lcoe_occ * annual_generation / crf

    @info "  Back-calculated OCC: $(round(occ_eur/1e6, digits=1)) million EUR"

    # Calculate IDC in EUR for base and delayed scenarios
    idc_base_eur = occ_eur * idc_factor(base_construction_time, wacc)
    idc_delayed_eur = occ_eur * idc_factor(delayed_construction_time, wacc)

    # Convert IDC back to EUR/MWh
    lcoe_idc_base = (idc_base_eur * crf) / annual_generation
    lcoe_idc_delayed = (idc_delayed_eur * crf) / annual_generation

    # Total capital LCOE for each scenario
    lcoe_capital_base = lcoe_occ + lcoe_idc_base
    lcoe_capital_delayed = lcoe_occ + lcoe_idc_delayed

    @info "  On-time ($(base_construction_time)yr): OCC=$(round(lcoe_occ, digits=1)), IDC=$(round(lcoe_idc_base, digits=1)), Total=$(round(lcoe_capital_base, digits=1)) EUR/MWh"
    @info "  Delayed ($(delayed_construction_time)yr): OCC=$(round(lcoe_occ, digits=1)), IDC=$(round(lcoe_idc_delayed, digits=1)), Total=$(round(lcoe_capital_delayed, digits=1)) EUR/MWh"

    # Create figure with vertically stacked layout for better readability in text
    fig = Figure(size=(900, 1100), backgroundcolor=:white)

    # Scenario 1: On-time (top panel)
    ax1 = Axis(fig[1,1],
              title = "(a) On-Time Construction: $(base_construction_time) years",
              titlesize = 18,
              ylabel = "Capital LCOE [EUR₂₀₂₅/MWh]",
              ylabelsize = 16,
              xticks = (1:3, ["OCC\n(Overnight Cost)", "IDC\n(Financing)", "Total\nCapital"]),
              xticklabelsize = 14,
              yticklabelsize = 14,
              xticklabelrotation = 0,
              xgridvisible = false,
              ygridvisible = true,
              ygridcolor = (:gray, 0.2),
              ygridstyle = :dash)

    barplot!(ax1, [1], [lcoe_occ], color=RGBf(0.0, 0.5, 0.5), width=0.6)  # Teal (matches LCOE breakdown)
    barplot!(ax1, [2], [lcoe_idc_base], offset=[lcoe_occ], color=RGBf(0.8, 0.4, 0.0), width=0.6)  # Darker orange
    barplot!(ax1, [3], [lcoe_capital_base], color=RGBf(0.0, 0.4, 0.4), width=0.6)  # Darker teal for total

    # Labels - increased font sizes
    text!(ax1, 1, lcoe_occ/2,
         text="$(round(lcoe_occ, digits=1))\n($(round(lcoe_occ/lcoe_capital_base*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=16, font=:bold)
    text!(ax1, 2, lcoe_occ + lcoe_idc_base/2,
         text="$(round(lcoe_idc_base, digits=1))\n($(round(lcoe_idc_base/lcoe_capital_base*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=16, font=:bold)
    text!(ax1, 3, lcoe_capital_base/2,
         text="$(round(lcoe_capital_base, digits=1))",
         align=(:center, :center), color=:white, fontsize=18, font=:bold)

    lines!(ax1, [1.3, 1.7], [lcoe_occ, lcoe_occ], color=:gray50, linestyle=:dash, linewidth=2)
    ylims!(ax1, 0, max(lcoe_capital_base, lcoe_capital_delayed) * 1.2)

    # Scenario 2: Delayed (bottom panel)
    ax2 = Axis(fig[2,1],
              title = "(b) Delayed Construction: $(delayed_construction_time) years (+$(delayed_construction_time - base_construction_time) years delay)",
              titlesize = 18,
              ylabel = "Capital LCOE [EUR₂₀₂₅/MWh]",
              ylabelsize = 16,
              xticks = (1:3, ["OCC\n(Overnight Cost)", "IDC\n(Financing)", "Total\nCapital"]),
              xticklabelsize = 14,
              yticklabelsize = 14,
              xticklabelrotation = 0,
              xgridvisible = false,
              ygridvisible = true,
              ygridcolor = (:gray, 0.2),
              ygridstyle = :dash)

    barplot!(ax2, [1], [lcoe_occ], color=RGBf(0.0, 0.5, 0.5), width=0.6)  # Teal (matches LCOE breakdown)
    barplot!(ax2, [2], [lcoe_idc_delayed], offset=[lcoe_occ], color=RGBf(0.8, 0.2, 0.2), width=0.6)  # Orangered (delayed scenario)
    barplot!(ax2, [3], [lcoe_capital_delayed], color=RGBf(0.6, 0.15, 0.15), width=0.6)  # Darker red for total

    # Labels - increased font sizes
    text!(ax2, 1, lcoe_occ/2,
         text="$(round(lcoe_occ, digits=1))\n($(round(lcoe_occ/lcoe_capital_delayed*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=16, font=:bold)
    text!(ax2, 2, lcoe_occ + lcoe_idc_delayed/2,
         text="$(round(lcoe_idc_delayed, digits=1))\n($(round(lcoe_idc_delayed/lcoe_capital_delayed*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=16, font=:bold)
    text!(ax2, 3, lcoe_capital_delayed/2,
         text="$(round(lcoe_capital_delayed, digits=1))",
         align=(:center, :center), color=:white, fontsize=18, font=:bold)

    lines!(ax2, [1.3, 1.7], [lcoe_occ, lcoe_occ], color=:gray50, linestyle=:dash, linewidth=2)
    ylims!(ax2, 0, max(lcoe_capital_base, lcoe_capital_delayed) * 1.2)

    # Overall title
    Label(fig[0, 1],
          "Construction Delay Impact on Capital Costs\n$(scale) Reactor @ $(round(wacc*100, digits=0))% WACC",
          fontsize=20, font=:bold)

    # Impact summary
    idc_increase = lcoe_idc_delayed - lcoe_idc_base
    total_increase = lcoe_capital_delayed - lcoe_capital_base
    pct_increase = (total_increase / lcoe_capital_base) * 100
    idc_per_year = idc_increase / (delayed_construction_time - base_construction_time)

    Label(fig[3, 1],
          "Impact: +$(round(idc_per_year, digits=1)) EUR/MWh per year of delay\nTotal: +$(round(total_increase, digits=1)) EUR/MWh (+$(round(pct_increase, digits=1))%)",
          fontsize=16, color=:darkred, font=:bold)

    # Legend (updated to match consistent color scheme)
    legend_elements = [
        PolyElement(color=RGBf(0.0, 0.5, 0.5), strokewidth=0),  # Teal
        PolyElement(color=RGBf(0.8, 0.4, 0.0), strokewidth=0),  # Darker orange
        PolyElement(color=RGBf(0.8, 0.2, 0.2), strokewidth=0)   # Orangered
    ]
    legend_labels = [
        "OCC (Overnight Cost)",
        "IDC - On-time",
        "IDC - Delayed"
    ]

    Legend(fig[4, 1],
           legend_elements,
           legend_labels,
           orientation=:horizontal,
           tellwidth=false,
           tellheight=true,
           framevisible=true,
           labelsize=14)

    return fig
end

"""
    plot_lcoe_breakdown_aggregate(pjs::Vector, scale::String, wacc_range::Vector,
                                   electricity_price_mean::Float64, opt_scaling::String,
                                   construction_time_range::Vector{Int})

Create aggregated LCOE breakdown visualization with both stacked bar and donut chart.
Shows the component breakdown (OCC, IDC, Fixed O&M, Variable O&M+Fuel) for a typical reactor of given scale.
"""
function plot_lcoe_breakdown_aggregate(pjs::Vector, scale::String, wacc_range::Vector,
                                       electricity_price_mean::Float64, opt_scaling::String,
                                       construction_time_range::Vector{Int})

    # FIXED: Use only BWRX-300 to match standalone donut and IDC sensitivity
    scale_reactors = filter(p -> p.name == "BWRX-300", pjs)
    if isempty(scale_reactors)
        @error "BWRX-300 reactor not found"
        return nothing
    end

    @info "Creating LCOE breakdown for BWRX-300 (using $opt_scaling scaling)"

    n = 1_000_000  # Match main MCS and standalone donut
    all_capital = Float64[]
    all_wacc = Float64[]
    all_construction_time = Int[]
    all_fixed_om = Float64[]
    all_variable_om_fuel = Float64[]
    all_total = Float64[]  # FIXED: Track total LCOE for correct median

    for p in scale_reactors
        rand_vars = gen_rand_vars(opt_scaling, n, wacc_range, electricity_price_mean, p;
                                  construction_time_range=construction_time_range)
        disc_res = mc_run(n, p, rand_vars)
        res = npv_lcoe(disc_res, decompose=true)

        # Collect capital cost and parameters for OCC/IDC decomposition
        append!(all_capital, res.lcoe_capital)
        append!(all_wacc, rand_vars.wacc)
        append!(all_construction_time, rand_vars.construction_time)
        append!(all_fixed_om, res.lcoe_fixed_om)
        append!(all_variable_om_fuel, res.lcoe_variable_om_fuel)
        # FIXED: Use authoritative total LCOE from npv_lcoe function
        append!(all_total, res.lcoe)
    end

    # Decompose capital into OCC and IDC (illustrative, using Wealer formula)
    occ_illustrative, idc_illustrative = decompose_capital_to_occ_idc(
        all_capital, all_wacc, all_construction_time
    )

    # Take medians
    median_occ = median(occ_illustrative)
    median_idc = median(idc_illustrative)
    median_fixed_om = median(all_fixed_om)
    median_variable = median(all_variable_om_fuel)
    median_total = median(all_total)  # FIXED: median of totals, not sum of medians

    @info "  Median LCOE: $(round(median_total, digits=1)) EUR2025/MWh"
    @info "    OCC: $(round(median_occ, digits=1)) ($(round(median_occ/median_total*100, digits=1))%)"
    @info "    IDC: $(round(median_idc, digits=1)) ($(round(median_idc/median_total*100, digits=1))%)"
    @info "    Fixed O&M: $(round(median_fixed_om, digits=1)) ($(round(median_fixed_om/median_total*100, digits=1))%)"
    @info "    Variable O&M+Fuel: $(round(median_variable, digits=1)) ($(round(median_variable/median_total*100, digits=1))%)"

    # Define colors for all 4 components (matching standalone donut)
    component_colors = Dict(
        "OCC" => RGBf(0.2, 0.4, 0.8),
        "IDC" => RGBf(0.9, 0.5, 0.2),
        "Fixed O&M" => RGBf(0.3, 0.7, 0.3),
        "Variable O&M + Fuel" => RGBf(0.9, 0.3, 0.3)
    )

    fig = Figure(size=(1400, 700), backgroundcolor=:white)

    ax1 = Axis(fig[1,1],
              title = "Stacked Bar Chart",
              ylabel = "LCOE [EUR2025/MWh]",
              xticks = ([1], ["BWRX-300"]),
              xgridvisible = false,
              ygridvisible = true,
              ygridcolor = (:gray, 0.2),
              ygridstyle = :dash)

    # Stack all 4 components
    barplot!(ax1, [1], [median_occ], color=component_colors["OCC"], width=0.6)
    barplot!(ax1, [1], [median_idc], offset=[median_occ],
             color=component_colors["IDC"], width=0.6)
    barplot!(ax1, [1], [median_fixed_om], offset=[median_occ + median_idc],
             color=component_colors["Fixed O&M"], width=0.6)
    barplot!(ax1, [1], [median_variable], offset=[median_occ + median_idc + median_fixed_om],
             color=component_colors["Variable O&M + Fuel"], width=0.6)

    # Add text labels to each component
    text!(ax1, 1, median_occ/2,
         text="OCC\n$(round(median_occ, digits=1))\n($(round(median_occ/median_total*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=12, font=:bold)

    text!(ax1, 1, median_occ + median_idc/2,
         text="IDC\n$(round(median_idc, digits=1))\n($(round(median_idc/median_total*100, digits=1))%)",
         align=(:center, :center), color=:black, fontsize=12, font=:bold)

    text!(ax1, 1, median_occ + median_idc + median_fixed_om/2,
         text="Fixed O&M\n$(round(median_fixed_om, digits=1))\n($(round(median_fixed_om/median_total*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=12, font=:bold)

    text!(ax1, 1, median_occ + median_idc + median_fixed_om + median_variable/2,
         text="Variable\n$(round(median_variable, digits=1))\n($(round(median_variable/median_total*100, digits=1))%)",
         align=(:center, :center), color=:white, fontsize=12, font=:bold)

    text!(ax1, 1, median_total + median_total*0.05,
         text="Total: $(round(median_total, digits=1))",
         align=(:center, :bottom), color=:black, fontsize=13, font=:bold)

    ax2 = Axis(fig[1,2], title = "Donut Chart", aspect = DataAspect())
    hidedecorations!(ax2)
    hidespines!(ax2)

    values = [median_occ, median_idc, median_fixed_om, median_variable]
    colors_donut = [component_colors["OCC"],
                    component_colors["IDC"],
                    component_colors["Fixed O&M"],
                    component_colors["Variable O&M + Fuel"]]
    labels = ["OCC", "IDC", "Fixed O&M", "Variable O&M+Fuel"]

    # FIXED: Use median_total (median of totals), not sum(values) (sum of medians)
    angles = [0.0]
    for v in values
        push!(angles, angles[end] + 2π * v / median_total)
    end

    for i in 1:length(values)
        θ_range = range(angles[i], angles[i+1], length=100)
        outer_radius = 1.0
        inner_radius = 0.5

        x_outer = outer_radius .* cos.(θ_range)
        y_outer = outer_radius .* sin.(θ_range)
        x_inner = inner_radius .* cos.(reverse(θ_range))
        y_inner = inner_radius .* sin.(reverse(θ_range))

        x_poly = vcat(x_outer, x_inner)
        y_poly = vcat(y_outer, y_inner)

        poly!(ax2, Point2f.(x_poly, y_poly), color=colors_donut[i], strokewidth=2, strokecolor=:white)

        mid_angle = (angles[i] + angles[i+1]) / 2
        label_radius = 0.75
        x_label = label_radius * cos(mid_angle)
        y_label = label_radius * sin(mid_angle)

        # FIXED: Calculate percentage against median_total, not sum(values)
        pct = round(values[i] / median_total * 100, digits=1)
        text!(ax2, x_label, y_label,
             text="$(labels[i])\n$(round(values[i], digits=1))\n($(pct)%)",
             align=(:center, :center), color=:white, fontsize=11, font=:bold)
    end

    xlims!(ax2, -1.3, 1.3)
    ylims!(ax2, -1.3, 1.3)

    Label(fig[0, :],
          "LCOE Component Breakdown: BWRX-300 Reactor ($opt_scaling scaling)",
          fontsize=20, font=:bold)

    legend_elements = [
        PolyElement(color=component_colors["OCC"], strokewidth=0),
        PolyElement(color=component_colors["IDC"], strokewidth=0),
        PolyElement(color=component_colors["Fixed O&M"], strokewidth=0),
        PolyElement(color=component_colors["Variable O&M + Fuel"], strokewidth=0)
    ]
    legend_labels = [
        "OCC (Overnight Construction Cost)",
        "IDC (Interest During Construction)",
        "Fixed O&M (Operations & Maintenance)",
        "Variable O&M + Fuel"
    ]

    Legend(fig[2, :],
           legend_elements,
           legend_labels,
           orientation=:horizontal,
           tellwidth=false,
           tellheight=true,
           framevisible=true)

    # Add caption explaining illustrative decomposition
    Label(fig[3, :],
          "Note: OCC and IDC shown as illustrative decomposition using Wealer formula: IDC = OCC × [(r/2)×T + (r²/6)×T²].\nSimulation uses DCF framework where IDC is implicit in present value calculations.",
          fontsize=10,
          tellwidth=false,
          justification=:left)

    return fig
end

"""
    plot_lcoe_donut_standalone(pjs::Vector, scale::String, wacc_range::Vector,
                                electricity_price_mean::Float64, opt_scaling::String,
                                construction_time_range::Vector{Int})

Create standalone donut chart for LCOE breakdown with vibrant colors and borders.
"""
function plot_lcoe_donut_standalone(pjs::Vector, scale::String, wacc_range::Vector,
                                    electricity_price_mean::Float64, opt_scaling::String,
                                    construction_time_range::Vector{Int})

    # Filter for BWRX-300 only
    scale_reactors = filter(p -> p.name == "BWRX-300", pjs)
    if isempty(scale_reactors)
        @error "BWRX-300 reactor not found"
        return nothing
    end

    @info "Creating standalone donut chart for BWRX-300"

    n = 1_000_000  # Match main MCS sample size for accuracy
    all_capital = Float64[]
    all_wacc = Float64[]
    all_construction_time = Int[]
    all_fixed_om = Float64[]
    all_variable_om_fuel = Float64[]
    all_total = Float64[]

    for p in scale_reactors
        rand_vars = gen_rand_vars(opt_scaling, n, wacc_range, electricity_price_mean, p;
                                  construction_time_range=construction_time_range)
        disc_res = mc_run(n, p, rand_vars)
        res = npv_lcoe(disc_res, decompose=true)

        # Collect capital cost and parameters for OCC/IDC decomposition
        append!(all_capital, res.lcoe_capital)
        append!(all_wacc, rand_vars.wacc)
        append!(all_construction_time, rand_vars.construction_time)
        append!(all_fixed_om, res.lcoe_fixed_om)
        append!(all_variable_om_fuel, res.lcoe_variable_om_fuel)
        append!(all_total, res.lcoe)
    end

    # Decompose capital into OCC and IDC (illustrative, using Wealer formula)
    occ_illustrative, idc_illustrative = decompose_capital_to_occ_idc(
        all_capital, all_wacc, all_construction_time
    )

    # Take medians
    median_occ = median(occ_illustrative)
    median_idc = median(idc_illustrative)
    median_fixed_om = median(all_fixed_om)
    median_variable = median(all_variable_om_fuel)
    median_total = median(all_total)

    @info "  Median LCOE: $(round(median_total, digits=1)) EUR2025/MWh"
    @info "    OCC (illustrative): $(round(median_occ, digits=1)) ($(round(median_occ/median_total*100, digits=1))%)"
    @info "    IDC (illustrative): $(round(median_idc, digits=1)) ($(round(median_idc/median_total*100, digits=1))%)"
    @info "    Fixed O&M: $(round(median_fixed_om, digits=1)) ($(round(median_fixed_om/median_total*100, digits=1))%)"
    @info "    Variable: $(round(median_variable, digits=1)) ($(round(median_variable/median_total*100, digits=1))%)"

    # More vibrant colors - 4 components
    component_colors = Dict(
        "OCC" => RGBf(0.0, 0.65, 0.65),                   # Brighter teal
        "IDC" => RGBf(1.0, 0.55, 0.0),                    # Vibrant orange
        "Fixed O&M" => RGBf(0.95, 0.85, 0.0),             # Bright yellow
        "Variable O&M + Fuel" => RGBf(1.0, 0.3, 0.3)      # Vibrant red
    )

    # Create figure - horizontal layout with legend on right
    fig = Figure(size=(1200, 600), backgroundcolor=:white)

    # Create main grid layout
    gl = fig[1, 1] = GridLayout()

    # Create donut chart on left
    ax = Axis(gl[1, 1], aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)

    # Four components: OCC, IDC (illustrative decomposition), Fixed O&M, Variable
    values = [median_occ, median_idc, median_fixed_om, median_variable]
    colors_donut = [component_colors["OCC"], component_colors["IDC"],
                    component_colors["Fixed O&M"], component_colors["Variable O&M + Fuel"]]

    # FIXED: Use median_total (median of totals), not sum(values) (sum of medians)
    # This ensures donut slices match the percentages shown in legend
    angles = [0.0]
    for v in values
        push!(angles, angles[end] + 2π * v / median_total)
    end

    # Draw donut segments with borders
    for i in 1:length(values)
        θ_range = range(angles[i], angles[i+1], length=100)
        outer_radius = 1.0
        inner_radius = 0.55  # Slightly thicker donut

        x_outer = outer_radius .* cos.(θ_range)
        y_outer = outer_radius .* sin.(θ_range)
        x_inner = inner_radius .* cos.(reverse(θ_range))
        y_inner = inner_radius .* sin.(reverse(θ_range))

        x_poly = vcat(x_outer, x_inner)
        y_poly = vcat(y_outer, y_inner)

        # Draw filled polygon with border
        poly!(ax, Point2f.(x_poly, y_poly),
             color=colors_donut[i],
             strokewidth=2,
             strokecolor=:white)
    end

    # Add center text showing total
    text!(ax, 0, 0,
         text="Total LCOE\n$(round(median_total, digits=1))\nEUR2025/MWh",
         align=(:center, :center),
         color=:black,
         fontsize=24,
         font=:bold)

    # Set axis limits
    xlims!(ax, -1.15, 1.15)
    ylims!(ax, -1.15, 1.15)

    # Legend on the right side with percentages
    legend_elements = [
        PolyElement(color=component_colors["OCC"], strokewidth=2, strokecolor=:black),
        PolyElement(color=component_colors["IDC"], strokewidth=2, strokecolor=:black),
        PolyElement(color=component_colors["Fixed O&M"], strokewidth=2, strokecolor=:black),
        PolyElement(color=component_colors["Variable O&M + Fuel"], strokewidth=2, strokecolor=:black)
    ]

    # Format labels with percentages
    pct_occ = round(median_occ / median_total * 100, digits=1)
    pct_idc = round(median_idc / median_total * 100, digits=1)
    pct_fixed = round(median_fixed_om / median_total * 100, digits=1)
    pct_variable = round(median_variable / median_total * 100, digits=1)

    legend_labels = [
        "OCC: $(round(median_occ, digits=1)) EUR/MWh ($(pct_occ)%)",
        "IDC: $(round(median_idc, digits=1)) EUR/MWh ($(pct_idc)%)",
        "Fixed O&M: $(round(median_fixed_om, digits=1)) EUR/MWh ($(pct_fixed)%)",
        "Variable O&M + Fuel: $(round(median_variable, digits=1)) EUR/MWh ($(pct_variable)%)"
    ]

    Legend(gl[1, 2],
           legend_elements,
           legend_labels,
           orientation=:vertical,
           tellwidth=true,
           tellheight=false,
           framevisible=true,
           labelsize=20,
           rowgap=15,
           patchsize=(25, 25),
           padding=(25, 25, 25, 25))

    # Add caption explaining illustrative decomposition - matching legend text size
    Label(fig[2, 1],
          "Note: OCC and IDC shown as illustrative decomposition using Wealer formula: IDC = OCC × [(r/2)×T + (r²/6)×T²].\nSimulation uses DCF framework where IDC is implicit in present value calculations.",
          fontsize=18,
          tellwidth=false,
          justification=:left)

    return fig
end

"""
    shapley_plot_single_scale(shapley_results, scale::String, pjs::Vector)

Create Shapley sensitivity heatmap for a single reactor scale.
Returns a figure showing only one scale (Micro, SMR, or Large).
"""
function shapley_plot_single_scale(shapley_results, scale::String, pjs::Vector)

    # Filter reactors by scale
    reactors_in_scale = [pj.name for pj in pjs if pj.scale == scale && pj.name in names(shapley_results)]

    if isempty(reactors_in_scale)
        @warn "No reactors found for scale: $scale"
        return nothing
    end

    # Extract variable names and map to display names
    var_display_names = Dict(
        "wacc" => "WACC",
        "construction_time" => "Construction Time",
        "capacity factor" => "Capacity Factor",
        "scaling" => "Investment OCC"
    )

    xticks_raw = shapley_results.var
    xticks = [get(var_display_names, var, var) for var in xticks_raw]

    # Extract data for this scale
    data_shapley = Matrix(shapley_results[:, reactors_in_scale])
    yticks_scale = reactors_in_scale

    # Sizing for readable cells with text values
    n_reactors = length(reactors_in_scale)
    n_vars = length(xticks)
    cell_width = 150  # Wide enough for text
    cell_height = 60
    type_bar_width = 1  # Minimal width to eliminate empty space

    fig_width = type_bar_width + (n_vars * cell_width) + 250  # Type bar + heatmap + colorbar
    fig_height = n_reactors * cell_height + 250

    # Create figure
    fig = Figure(size = (fig_width, fig_height), backgroundcolor=:white)
    colormap_shapley = :deep  # Same as original uncertainty report

    # Main grid layout
    gl = fig[1, 1] = GridLayout()

    # Reactor type indicators
    type_indicator_colors_map = Dict(
        "PWR" => RGBf(0.9, 0.5, 0.1),
        "BWR" => RGBf(0.9, 0.5, 0.1),
        "HTR" => RGBf(0.9, 0.8, 0.1),
        "SFR" => RGBf(0.1, 0.7, 0.7),
        "MSR" => RGBf(0.6, 0.6, 0.6)
    )

    reactor_types_in_scale = []
    for reactor_name in reactors_in_scale
        pj_idx = findfirst(p -> p.name == reactor_name, pjs)
        if !isnothing(pj_idx)
            push!(reactor_types_in_scale, pjs[pj_idx].type)
        else
            push!(reactor_types_in_scale, "Unknown")
        end
    end

    # Shapley heatmap (no separate type column - circles will be on heatmap)
    ax_shapley = Axis(gl[1, 1],
                      xticks = (1:length(xticks), xticks),
                      yticks = (1:length(yticks_scale), yticks_scale),
                      title = "Shapley Sensitivity Indices",
                      titlesize = 18,
                      yaxisposition = :right,
                      yticklabelsize = 16,
                      xticklabelsize = 14)
    ax_shapley.xticklabelrotation = π / 3
    ax_shapley.xticklabelalign = (:right, :center)

    # Extend x-axis limits to show reactor type indicator circles
    xlims!(ax_shapley, 0, length(xticks) + 0.5)

    # Draw heatmap
    hm_shapley = heatmap!(ax_shapley, 1:length(xticks), 1:length(yticks_scale), data_shapley,
                          colormap = colormap_shapley,
                          colorrange = (0, 1))

    # Add text annotations showing values
    for i in 1:n_vars, j in 1:n_reactors
        val = data_shapley[i, j]
        text!(ax_shapley, i, j, text = string(round(val, digits=2)),
              align = (:center, :center),
              fontsize = 14,
              color = val > 0.5 ? :white : :black)
    end

    # Draw reactor type indicators as circles on the left side of heatmap
    for (i, reactor_type) in enumerate(reactor_types_in_scale)
        type_color = get(type_indicator_colors_map, reactor_type, RGBf(0.5, 0.5, 0.5))
        scatter!(ax_shapley, [0.3], [i],
                color=type_color,
                strokecolor=:black,
                strokewidth=1.5,
                markersize=15)
    end

    # Colorbar
    Colorbar(gl[1, 2], hm_shapley, label = "Shapley Value", labelsize = 16, ticklabelsize = 14)

    # Legend for reactor types (bottom)
    unique_types = unique(reactor_types_in_scale)
    legend_elements = [PolyElement(color=type_indicator_colors_map[t], strokewidth=1, strokecolor=:black)
                       for t in unique_types if haskey(type_indicator_colors_map, t)]
    legend_labels = unique_types[findall(t -> haskey(type_indicator_colors_map, t), unique_types)]

    Legend(fig[2, 1], legend_elements, legend_labels,
           "Reactor Type",
           orientation = :horizontal,
           framevisible = false,
           tellwidth = false,
           tellheight = true,
           labelsize = 14,
           titlesize = 16)

    return fig
end

"""
    si_plot_single_scale(si_results, scale::String, pjs::Vector)

Create Sobol sensitivity heatmap for a single reactor scale.
Returns a figure showing only one scale (Micro, SMR, or Large).
"""
function si_plot_single_scale(si_results, scale::String, pjs::Vector)

    # Filter reactors by scale
    reactors_in_scale = [pj.name for pj in pjs if pj.scale == scale && pj.name in names(si_results)]

    if isempty(reactors_in_scale)
        @warn "No reactors found for scale: $scale"
        return nothing
    end

    # Extract variable names
    var_display_names = Dict(
        "wacc" => "WACC",
        "construction_time" => "Construction Time",
        "capacity factor" => "Capacity Factor",
        "scaling" => "Investment OCC"
    )

    xticks_raw = unique(si_results.var)
    xticks = [get(var_display_names, var, var) for var in xticks_raw]

    # Extract data for this scale
    si_s = filter(:si => index -> index == "S", si_results)
    si_st = filter(:si => index -> index == "ST", si_results)

    data_s = Matrix(si_s[:, reactors_in_scale])
    data_st = Matrix(si_st[:, reactors_in_scale])

    yticks_scale = reactors_in_scale

    # Sizing: side-by-side layout needs width for two heatmaps + colorbar + legend
    n_reactors = length(reactors_in_scale)
    n_vars = length(xticks)
    cell_width = 120  # Smaller cells for better text-to-box ratio
    cell_height = 50   # Smaller rows for better text-to-box ratio
    type_bar_width = 1  # Minimal width to eliminate empty space

    fig_width = type_bar_width + 2 * (n_vars * cell_width) + 350  # Two heatmaps + colorbar + margins
    fig_height = n_reactors * cell_height + 280  # Reactor rows + title + legend

    # Create figure with side-by-side layout
    fig = Figure(size = (fig_width, fig_height), backgroundcolor=:white)
    colormap_si = :deep  # Same as original uncertainty report

    # Main grid layout
    gl = fig[1, 1] = GridLayout()

    # Reactor type indicators (shared on left)
    type_indicator_colors_map = Dict(
        "PWR" => RGBf(0.9, 0.5, 0.1),
        "BWR" => RGBf(0.9, 0.5, 0.1),
        "HTR" => RGBf(0.9, 0.8, 0.1),
        "SFR" => RGBf(0.1, 0.7, 0.7),
        "MSR" => RGBf(0.6, 0.6, 0.6)
    )

    reactor_types_in_scale = []
    for reactor_name in reactors_in_scale
        pj_idx = findfirst(p -> p.name == reactor_name, pjs)
        if !isnothing(pj_idx)
            push!(reactor_types_in_scale, pjs[pj_idx].type)
        else
            push!(reactor_types_in_scale, "Unknown")
        end
    end

    # First-order indices (S) - column 1
    ax_s = Axis(gl[1, 1],
                xticks = (1:length(xticks), xticks),
                yticks = (1:length(yticks_scale), fill("", n_reactors)),  # Empty labels (shown on left)
                title = "First-order Effect",
                titlesize = 18,
                xticklabelsize = 14)
    ax_s.xticklabelrotation = π / 3
    ax_s.xticklabelalign = (:right, :center)
    ax_s.yticklabelsvisible = false

    # Extend x-axis limits to show reactor type indicator circles
    xlims!(ax_s, 0, length(xticks) + 0.5)

    # Draw S heatmap
    hm_s = heatmap!(ax_s, 1:length(xticks), 1:length(yticks_scale), data_s,
                    colormap = colormap_si, colorrange = (0, 1))

    # Add text annotations showing values
    for i in 1:n_vars, j in 1:n_reactors
        val = data_s[i, j]
        text!(ax_s, i, j, text = string(round(val, digits=2)),
              align = (:center, :center),
              fontsize = 16,
              font = :bold,
              color = val > 0.5 ? :white : :black)
    end

    # Draw reactor type indicators as circles on the left side of S heatmap
    for (i, reactor_type) in enumerate(reactor_types_in_scale)
        type_color = get(type_indicator_colors_map, reactor_type, RGBf(0.5, 0.5, 0.5))
        scatter!(ax_s, [0.3], [i],
                color=type_color,
                strokecolor=:black,
                strokewidth=1.5,
                markersize=15)
    end

    # Total-order indices (ST) - column 2
    ax_st = Axis(gl[1, 2],
                 xticks = (1:length(xticks), xticks),
                 yticks = (1:length(yticks_scale), yticks_scale),  # Show labels on right
                 title = "Total-order Effect",
                 titlesize = 18,
                 xticklabelsize = 14,
                 yaxisposition = :right)
    ax_st.xticklabelrotation = π / 3
    ax_st.xticklabelalign = (:right, :center)
    ax_st.yticklabelsize = 16

    # Draw ST heatmap
    hm_st = heatmap!(ax_st, 1:length(xticks), 1:length(yticks_scale), data_st,
                     colormap = colormap_si, colorrange = (0, 1))

    # Add text annotations showing values
    for i in 1:n_vars, j in 1:n_reactors
        val = data_st[i, j]
        text!(ax_st, i, j, text = string(round(val, digits=2)),
              align = (:center, :center),
              fontsize = 16,
              font = :bold,
              color = val > 0.5 ? :white : :black)
    end

    # Shared colorbar
    Colorbar(gl[1, 3], hm_s, label = "Sensitivity Index", labelsize = 16, ticklabelsize = 14)

    # Legend for reactor types (bottom)
    unique_types = unique(reactor_types_in_scale)
    legend_elements = [PolyElement(color=type_indicator_colors_map[t], strokewidth=1, strokecolor=:black)
                       for t in unique_types if haskey(type_indicator_colors_map, t)]
    legend_labels = unique_types[findall(t -> haskey(type_indicator_colors_map, t), unique_types)]

    Legend(fig[2, 1], legend_elements, legend_labels,
           "Reactor Type",
           orientation = :horizontal,
           framevisible = false,
           tellwidth = false,
           tellheight = true,
           labelsize = 14,
           titlesize = 16)

    return fig
end
"""
    create_idc_sensitivity_table(scale_pjs, wacc_values, construction_times, opt_scaling)

Create a sensitivity table showing IDC component of LCOE as a function of construction time
and WACC, using a fixed median OCC value. Similar to Table 8.1 from NEA reference paper.

# Arguments
- `scale_pjs::Vector{project}`: Reactors of a specific scale (e.g., all SMRs)
- `wacc_values::Vector{Float64}`: WACC rates to test (e.g., [0.04, 0.07, 0.10])
- `construction_times::Vector{Int}`: Construction time scenarios (e.g., [3, 5, 7])
- `opt_scaling::String`: Scaling method for OCC calculation

# Returns
- DataFrame with IDC values (EUR/MWh) organized by construction time (rows) and WACC (columns)
- Figure with formatted table visualization
"""
function create_idc_sensitivity_table(scale_pjs::Vector,
                                     wacc_values::Vector{Float64},
                                     construction_times::Vector{Int},
                                     opt_scaling::String)

    scale_name = scale_pjs[1].scale
    n_reactors = length(scale_pjs)

    @info "Creating IDC sensitivity table for $(scale_name) reactors (n=$(n_reactors))"
    @info "  Construction times: $(construction_times) years"
    @info "  WACC values: $(wacc_values .* 100)%"

    # Calculate median OCC across all reactors (using manufacturer values with scaling)
    occ_values = Float64[]
    capacity_values = Float64[]
    loadfactor_values = Float64[]

    for pj in scale_pjs
        # Apply scaling to get median OCC
        scaling_factor = mean([0.4, 0.7])  # Use midpoint of scaling range
        scaled_investment = pj.investment * scaling_factor
        occ = scaled_investment * pj.plant_capacity  # EUR
        push!(occ_values, occ)
        push!(capacity_values, pj.plant_capacity)
        push!(loadfactor_values, mean(pj.loadfactor))
    end

    # Use median OCC for all calculations
    median_occ = median(occ_values)
    median_capacity = median(capacity_values)
    median_loadfactor = median(loadfactor_values)
    median_lifetime = scale_pjs[1].time[2]  # Operating years (should be same for all)

    @info "  Using median OCC: $(round(median_occ/1e6, digits=1)) M EUR"
    @info "  Median capacity: $(round(median_capacity, digits=0)) MWe"
    @info "  Median load factor: $(round(median_loadfactor*100, digits=1))%"

    # Annual generation [MWh/year]
    annual_generation = median_capacity * median_loadfactor * 8760

    # Initialize results DataFrame
    results = DataFrame(
        construction_time = Int[],
        wacc_4pct = Float64[],
        wacc_7pct = Float64[],
        wacc_10pct = Float64[]
    )

    # Calculate IDC for each combination
    for t_con in construction_times
        idc_by_wacc = Float64[]

        for wacc in wacc_values
            # Calculate IDC using explicit formula
            # IDC = OCC × [(r/2)×T_con + (r²/6)×T_con²]
            idc_factor = (wacc/2) * t_con + (wacc^2/6) * t_con^2
            idc_total = median_occ * idc_factor  # EUR

            # Calculate capital recovery factor
            crf = wacc * (1 + wacc)^median_lifetime / ((1 + wacc)^median_lifetime - 1)

            # Convert IDC to LCOE component [EUR/MWh]
            idc_lcoe = (idc_total * crf) / annual_generation

            push!(idc_by_wacc, idc_lcoe)

            @info "    t=$(t_con)yr, WACC=$(round(wacc*100, digits=0))%: IDC = $(round(idc_lcoe, digits=1)) EUR/MWh"
        end

        # Store results (assuming wacc_values = [0.04, 0.07, 0.10])
        push!(results, (t_con, idc_by_wacc[1], idc_by_wacc[2], idc_by_wacc[3]))
    end

    # Create formatted table visualization
    fig = Figure(resolution = (800, 500))

    # Add title
    Label(fig[1, 1],
          "IDC Component of LCOE - $(scale_name) Reactors",
          fontsize = 16, font = :bold, tellwidth = false)

    # Subtitle with key parameters
    Label(fig[2, 1],
          "Fixed OCC = $(round(median_occ/1e6, digits=1)) M EUR | Capacity = $(round(median_capacity, digits=0)) MWe | Load Factor = $(round(median_loadfactor*100, digits=1))%",
          fontsize = 11, tellwidth = false)

    # Create axis for table
    ax = Axis(fig[3, 1],
              xlabel = "",
              ylabel = "",
              aspect = DataAspect())

    hidedecorations!(ax)
    hidespines!(ax)

    # Table dimensions
    n_rows = length(construction_times) + 2  # +2 for header rows
    n_cols = 4  # Construction time + 3 WACC columns

    # Draw table grid (bold outer border)
    for i in 0:n_rows
        lw = (i == 0 || i == n_rows || i == 2) ? 2 : 1
        lines!(ax, [0, n_cols], [i, i], color = :black, linewidth = lw)
    end
    for j in 0:n_cols
        lw = (j == 0 || j == n_cols) ? 2 : 1
        lines!(ax, [j, j], [0, n_rows], color = :black, linewidth = lw)
    end

    # Add header - row 1 (main title)
    text!(ax, 2.0, n_rows - 0.5,
          text = "IDC Component of LCOE (EUR/MWh)",
          align = (:center, :center),
          fontsize = 12,
          font = :bold)

    # Add header - row 2 (column labels with shading)
    poly!(ax, [Point2f(0, n_rows-2), Point2f(n_cols, n_rows-2),
               Point2f(n_cols, n_rows-1), Point2f(0, n_rows-1)],
          color = RGBf(0.9, 0.9, 0.9))

    headers = ["Construction\nTime (years)", "WACC = 4%", "WACC = 7%", "WACC = 10%"]
    for (j, header) in enumerate(headers)
        text!(ax, j - 0.5, n_rows - 1.5,
              text = header,
              align = (:center, :center),
              fontsize = 11,
              font = :bold)
    end

    # Add data rows with alternating shading
    for (i, row) in enumerate(eachrow(results))
        row_idx = n_rows - i - 1

        # Alternating row shading
        if i % 2 == 0
            poly!(ax, [Point2f(0, row_idx-1), Point2f(n_cols, row_idx-1),
                       Point2f(n_cols, row_idx), Point2f(0, row_idx)],
                  color = RGBf(0.95, 0.95, 0.95))
        end

        # Construction time (bold)
        text!(ax, 0.5, row_idx - 0.5,
              text = string(row.construction_time),
              align = (:center, :center),
              fontsize = 11,
              font = :bold)

        # IDC LCOE values
        text!(ax, 1.5, row_idx - 0.5,
              text = string(round(row.wacc_4pct, digits=1)),
              align = (:center, :center),
              fontsize = 11)

        text!(ax, 2.5, row_idx - 0.5,
              text = string(round(row.wacc_7pct, digits=1)),
              align = (:center, :center),
              fontsize = 11)

        text!(ax, 3.5, row_idx - 0.5,
              text = string(round(row.wacc_10pct, digits=1)),
              align = (:center, :center),
              fontsize = 11)
    end

    # Add note
    Label(fig[4, 1],
          "Note: Values show the Interest During Construction (IDC) component of LCOE using median OCC from $(n_reactors) $(scale_name) reactors.\nIDC calculated using formula: IDC = OCC × [(r/2)×T + (r²/6)×T²], where r=WACC and T=construction time.\nLifetime = $(median_lifetime) years. Scaling method: $(opt_scaling).",
          fontsize = 9,
          tellwidth = false,
          justification = :left)

    xlims!(ax, 0, n_cols)
    ylims!(ax, 0, n_rows)

    return results, fig
end

"""
    plot_wacc_sensitivity_with_ci(pjs, wacc_range, opt_scaling, electricity_price_mean, construction_time_ranges)

Create WACC sensitivity plot showing median LCOE with confidence intervals,
similar to Figure 5.1 from NEA reference paper.

# Arguments
- `pjs::Vector{project}`: Reactor projects to analyze
- `wacc_range::Vector{Float64}`: Range of WACC values (e.g., 0.03:0.01:0.15)
- `opt_scaling::String`: Scaling method
- `electricity_price_mean::Float64`: Mean electricity price
- `construction_time_ranges::Dict`: Construction time ranges by scale

# Returns
- Figure with median LCOE lines and confidence interval bands for each reactor
"""
function plot_wacc_sensitivity_with_ci(pjs::Vector,
                                       wacc_range::AbstractVector{Float64},
                                       opt_scaling::String,
                                       electricity_price_mean::Float64,
                                       construction_time_ranges::Dict)

    @info "Creating WACC sensitivity plot with confidence intervals"
    @info "  WACC range: $(round(minimum(wacc_range)*100, digits=0))% - $(round(maximum(wacc_range)*100, digits=0))%"
    @info "  Reactors: $(length(pjs))"

    # Number of simulations per WACC value
    n = 10000

    # Store results for each reactor
    results_by_reactor = Dict{String, Dict{String, Vector{Float64}}}()

    # Run simulations for each reactor at each WACC value
    for pj in pjs
        @info "  Analyzing $(pj.name)..."

        median_lcoe = Float64[]
        p10_lcoe = Float64[]
        p90_lcoe = Float64[]

        for wacc in wacc_range
            # Generate random variables with small epsilon to avoid triangular dist errors
            # When min==max, add small range around the value
            epsilon = 0.001
            wacc_min = max(0.01, wacc - epsilon)
            wacc_max = min(0.20, wacc + epsilon)

            rand_vars = gen_rand_vars(opt_scaling, n, [wacc_min, wacc_max], electricity_price_mean, pj;
                                     construction_time_range=construction_time_ranges[pj.scale])

            # Run Monte Carlo simulation
            disc_res = mc_run(n, pj, rand_vars)

            # Calculate LCOE
            res = npv_lcoe(disc_res, decompose=false)

            # Calculate statistics
            push!(median_lcoe, median(res.lcoe))
            push!(p10_lcoe, quantile(res.lcoe, 0.10))
            push!(p90_lcoe, quantile(res.lcoe, 0.90))
        end

        results_by_reactor[pj.name] = Dict(
            "median" => median_lcoe,
            "p10" => p10_lcoe,
            "p90" => p90_lcoe
        )
    end

    # Create figure
    fig = Figure(resolution = (1000, 600))

    ax = Axis(fig[1, 1],
              xlabel = "Discount Rate (WACC)",
              ylabel = "LCOE (EUR/MWh)",
              title = "LCOE Sensitivity to Discount Rate with Confidence Intervals",
              titlesize = 14)

    # Color scheme by reactor type
    type_colors = Dict(
        "LWR" => RGBf(0.2, 0.4, 0.8),
        "SFR" => RGBf(0.8, 0.3, 0.3),
        "MSR" => RGBf(0.3, 0.7, 0.3),
        "HTGR" => RGBf(0.9, 0.6, 0.2),
        "micro" => RGBf(0.6, 0.4, 0.8)
    )

    wacc_pct = wacc_range .* 100  # Convert to percentage for x-axis

    # Plot each reactor
    for pj in pjs
        color = get(type_colors, pj.type, RGBf(0.5, 0.5, 0.5))

        median_vals = results_by_reactor[pj.name]["median"]
        p10_vals = results_by_reactor[pj.name]["p10"]
        p90_vals = results_by_reactor[pj.name]["p90"]

        # Plot confidence band (10th-90th percentile)
        band!(ax, wacc_pct, p10_vals, p90_vals,
              color = (color, 0.2),
              label = nothing)

        # Plot median line
        lines!(ax, wacc_pct, median_vals,
               color = color,
               linewidth = 2,
               label = pj.name)
    end

    # Add legend
    Legend(fig[1, 2], ax, "Reactor", framevisible = true)

    # Format x-axis as percentage
    ax.xticks = (0:5:20, ["$(i)%" for i in 0:5:20])

    return fig
end

"""
    calculate_vendor_baseline_lcoe_simple(pj, wacc_optimistic, construction_time_optimistic)

Calculate vendor's claimed LCOE using their stated OCC and OPTIMISTIC operating parameters.
This represents the LCOE at the manufacturer's claimed costs (N=1, no learning applied).

Uses optimistic assumptions to represent vendor's claimed competitiveness:
- Vendor-specified OCC (from manufacturer data)
- Fixed construction time at 3 years OR planned time (no delays)
- WACC at 7% (favorable financing conditions)
- Optimistic load factor (upper end of range)

# Arguments
- `pj::project`: Project object with reactor specifications
- `wacc_optimistic::Float64`: Optimistic discount rate (default 0.07 = 7%)
- `construction_time_optimistic::Int`: Optimistic construction time (default 3 years)

# Returns
- `vendor_baseline_lcoe::Float64`: LCOE at manufacturer's claimed costs
"""
function calculate_vendor_baseline_lcoe_simple(pj::project,
                                               wacc_optimistic::Float64=0.07,
                                               construction_time_optimistic::Int=3)
    # Use manufacturer OCC (no learning factor applied)
    occ = pj.investment * pj.plant_capacity  # EUR

    # Operating parameters (optimistic - use upper end of load factor range)
    loadfactor = maximum(pj.loadfactor)  # Optimistic load factor
    annual_generation = pj.plant_capacity * loadfactor * 8760  # MWh/year
    lifetime = pj.time[2]

    # Capital Recovery Factor
    crf = wacc_optimistic * (1 + wacc_optimistic)^lifetime / ((1 + wacc_optimistic)^lifetime - 1)

    # LCOE components
    # Capital (including IDC using explicit formula for vendor baseline)
    idc_factor = (wacc_optimistic/2) * construction_time_optimistic +
                 (wacc_optimistic^2/6) * construction_time_optimistic^2
    tcc = occ * (1 + idc_factor)
    lcoe_capital = (tcc * crf) / annual_generation

    # O&M
    fixed_om_annual = pj.plant_capacity * pj.operating_cost[1]
    variable_om_fuel = pj.operating_cost[2] + pj.operating_cost[3]

    lcoe_fixed_om = fixed_om_annual / annual_generation
    lcoe_variable = variable_om_fuel

    # Total
    vendor_baseline_lcoe = lcoe_capital + lcoe_fixed_om + lcoe_variable

    return vendor_baseline_lcoe
end

"""
    learning_curve_single_reactor(pj, wacc_range, electricity_price_mean, opt_scaling; kwargs...)

Create learning curve plot for a single reactor with vendor baseline.

# Arguments
- `pj::project`: Project object for reactor
- `wacc_range::Vector`: WACC range [min, max]
- `electricity_price_mean::Float64`: Mean electricity price
- `opt_scaling::String`: Scaling method identifier
- `N_max::Int=30`: Maximum number of cumulative units
- `learning_rates::Vector{Float64}=[0.05, 0.10, 0.15]`: Learning rates to plot
- `construction_time_range::Union{Nothing,Vector}=nothing`: Construction time range

# Returns
- `fig`: Figure with learning curves
- `vendor_baseline::Float64`: Vendor baseline LCOE
- `crossing_idx::Union{Int, Nothing}`: Index where base case crosses vendor baseline
"""
function learning_curve_single_reactor(pj::project, wacc_range::Vector,
                                      electricity_price_mean::Float64,
                                      opt_scaling::String;
                                      N_max::Int=30,
                                      learning_rates::Vector{Float64}=[0.05, 0.10, 0.15],
                                      construction_time_range::Union{Nothing,Vector}=nothing)

    @info "Creating learning curve for $(pj.name)"

    n_sim = 10000  # MC samples per learning scenario
    N_values = 1:N_max

    # Calculate vendor baseline using OPTIMISTIC assumptions:
    # - WACC = 7% (favorable financing)
    # - Construction time = 3 years (or planned time, no delays)
    # - Optimistic load factor (upper end of range)
    vendor_baseline = calculate_vendor_baseline_lcoe_simple(pj, 0.07, 3)

    @info "  Vendor baseline LCOE (optimistic): $(round(vendor_baseline, digits=1)) EUR/MWh"
    @info "    (WACC=7%, construction=3yrs, load factor=$(round(maximum(pj.loadfactor), digits=2)))"

    # Determine learning type based on scale
    learning_type = pj.scale == "Large" ? "deployment" : "factory"

    # Store results for each learning rate
    results = Dict{Float64, Vector{Float64}}()

    for LR in learning_rates
        lcoe_trajectory = Float64[]

        for N in N_values
            # Generate random variables with learning
            rand_vars = gen_rand_vars(opt_scaling, n_sim, wacc_range, electricity_price_mean, pj;
                                     apply_learning=true,
                                     N_unit=N,
                                     LR=LR,
                                     kappa=1.0,
                                     apply_soak_discount=false,
                                     construction_time_range=construction_time_range,
                                     quiet=true)

            disc_res = mc_run(n_sim, pj, rand_vars)
            res = npv_lcoe(disc_res, decompose=false)

            push!(lcoe_trajectory, median(res.lcoe))
        end

        results[LR] = lcoe_trajectory
    end

    # Create plot
    fig = Figure(size=(600, 500))

    ax = Axis(fig[1, 1],
              xlabel = "Cumulative Units (N)",
              ylabel = "LCOE [EUR2025/MWh]",
              title = pj.name,
              titlesize = 16,
              titlefont = :bold)

    # Color scheme
    lr_colors = Dict(
        0.05 => RGBf(0.7, 0.7, 0.7),     # Light gray
        0.10 => RGBf(0.2, 0.4, 0.8),     # Blue (base case)
        0.15 => RGBf(0.3, 0.7, 0.3)      # Green
    )

    lr_styles = Dict(
        0.05 => :dash,
        0.10 => :solid,
        0.15 => :dot
    )

    # Plot learning curves
    for LR in learning_rates
        lines!(ax, N_values, results[LR],
               color = lr_colors[LR],
               linestyle = lr_styles[LR],
               linewidth = 2.5,
               label = "LR $(Int(LR*100))%$(LR == 0.10 ? " (Base)" : "")")
    end

    # Add vendor baseline as horizontal line
    hlines!(ax, [vendor_baseline],
            color = :red,
            linestyle = :dashdot,
            linewidth = 2.5,
            label = "Vendor Baseline")

    # Find crossing points with vendor baseline for base case (LR=10%)
    base_trajectory = results[0.10]
    crossing_idx = findfirst(x -> x <= vendor_baseline, base_trajectory)

    if !isnothing(crossing_idx)
        crossing_n = N_values[crossing_idx]
        crossing_lcoe = base_trajectory[crossing_idx]

        # Add marker at crossing point
        scatter!(ax, [crossing_n], [crossing_lcoe],
                color = :red,
                markersize = 12,
                marker = :circle)

        # Add annotation
        text!(ax, crossing_n, crossing_lcoe,
              text = "  N=$(crossing_n)\n  $(round(crossing_lcoe, digits=1)) EUR/MWh",
              align = (:left, :bottom),
              fontsize = 11,
              color = :red,
              font = :bold)

        @info "  Crosses vendor baseline at N=$(crossing_n) units"
    else
        @info "  Does not reach vendor baseline within N=$(N_max) units"
    end

    # Legend
    Legend(fig[1, 1],
           ax,
           "",
           tellheight = false,
           tellwidth = false,
           halign = :right,
           valign = :top,
           margin = (10, 10, 10, 10),
           framevisible = true)

    # Add reactor type badge
    type_colors = Dict(
        "PWR" => RGBf(0.2, 0.4, 0.8),
        "BWR" => RGBf(0.3, 0.6, 0.9),
        "HTR" => RGBf(0.9, 0.6, 0.2),
        "SFR" => RGBf(0.8, 0.3, 0.3)
    )

    type_color = get(type_colors, pj.type, RGBf(0.5, 0.5, 0.5))

    Box(fig[1, 1],
        color = (type_color, 0.8),
        cornerradius = 5,
        strokewidth = 0,
        width = 50,
        height = 25,
        halign = :left,
        valign = :bottom,
        padding = (5, 5, 5, 5))

    Label(fig[1, 1],
          pj.type,
          fontsize = 11,
          font = :bold,
          color = :white,
          halign = :left,
          valign = :bottom,
          padding = (10, 10, 10, 10))

    return fig, vendor_baseline, crossing_idx
end

"""
    learning_curves_grid_appendix(pjs, wacc_range, electricity_price_mean, opt_scaling, construction_time_ranges)

Create 2-column × 5-row grid of learning curves for appendix.

# Arguments
- `pjs::Vector`: Vector of project objects
- `wacc_range::Vector`: WACC range [min, max]
- `electricity_price_mean::Float64`: Mean electricity price
- `opt_scaling::String`: Scaling method identifier
- `construction_time_ranges::Dict`: Construction time ranges by scale

# Returns
- `fig`: Figure with grid layout
- `crossing_data::DataFrame`: DataFrame with crossing point data
"""
function learning_curves_grid_appendix(pjs::Vector, wacc_range::Vector,
                                      electricity_price_mean::Float64,
                                      opt_scaling::String,
                                      construction_time_ranges::Dict)

    @info "Creating learning curves grid for appendix (2×5 layout)"

    # Create large figure for appendix
    fig = Figure(size=(1400, 3000))  # Tall figure

    n_reactors = length(pjs)
    n_cols = 2
    n_rows = ceil(Int, n_reactors / n_cols)

    crossing_data = DataFrame(
        Reactor = String[],
        VendorBaseline = Float64[],
        CrossingN = Union{Int, Missing}[],
        CrossingLCOE = Union{Float64, Missing}[]
    )

    idx = 1
    for row in 1:n_rows
        for col in 1:n_cols
            if idx > n_reactors
                break
            end

            pj = pjs[idx]

            @info "  Processing $(pj.name) ($(idx)/$(n_reactors))"

            # Calculate vendor baseline using OPTIMISTIC assumptions
            vendor_baseline = calculate_vendor_baseline_lcoe_simple(pj, 0.07, 3)

            # Run learning simulation
            n_sim = 10000
            N_max = 30
            learning_rates = [0.05, 0.10, 0.15]
            N_values = 1:N_max

            learning_type = pj.scale == "Large" ? "deployment" : "factory"

            results = Dict{Float64, Vector{Float64}}()

            for LR in learning_rates
                lcoe_trajectory = Float64[]

                for N in N_values
                    rand_vars = gen_rand_vars(opt_scaling, n_sim, wacc_range,
                                             electricity_price_mean, pj;
                                             apply_learning=true,
                                             N_unit=N,
                                             LR=LR,
                                             kappa=1.0,
                                             apply_soak_discount=false,
                                             construction_time_range=construction_time_ranges[pj.scale],
                                             quiet=true)

                    disc_res = mc_run(n_sim, pj, rand_vars)
                    res = npv_lcoe(disc_res, decompose=false)

                    push!(lcoe_trajectory, median(res.lcoe))
                end

                results[LR] = lcoe_trajectory
            end

            # Create subplot
            ax = Axis(fig[row, col],
                     xlabel = "Cumulative Units (N)",
                     ylabel = "LCOE [EUR2025/MWh]",
                     title = pj.name,
                     titlesize = 14,
                     titlefont = :bold)

            # Colors and styles
            lr_colors = Dict(
                0.05 => RGBf(0.7, 0.7, 0.7),
                0.10 => RGBf(0.2, 0.4, 0.8),
                0.15 => RGBf(0.3, 0.7, 0.3)
            )

            lr_styles = Dict(
                0.05 => :dash,
                0.10 => :solid,
                0.15 => :dot
            )

            # Plot curves
            for LR in learning_rates
                lines!(ax, N_values, results[LR],
                       color = lr_colors[LR],
                       linestyle = lr_styles[LR],
                       linewidth = 2)
            end

            # Vendor baseline
            hlines!(ax, [vendor_baseline],
                    color = :red,
                    linestyle = :dashdot,
                    linewidth = 2)

            # Find crossing
            base_trajectory = results[0.10]
            crossing_idx = findfirst(x -> x <= vendor_baseline, base_trajectory)

            if !isnothing(crossing_idx)
                crossing_n = N_values[crossing_idx]
                crossing_lcoe = base_trajectory[crossing_idx]

                scatter!(ax, [crossing_n], [crossing_lcoe],
                        color = :red,
                        markersize = 10,
                        marker = :circle)

                text!(ax, crossing_n, crossing_lcoe,
                      text = "  N=$(crossing_n)",
                      align = (:left, :bottom),
                      fontsize = 9,
                      color = :red,
                      font = :bold)

                push!(crossing_data, (pj.name, vendor_baseline, crossing_n, crossing_lcoe))
            else
                push!(crossing_data, (pj.name, vendor_baseline, missing, missing))
            end

            # Type badge
            type_colors = Dict(
                "PWR" => RGBf(0.2, 0.4, 0.8),
                "BWR" => RGBf(0.3, 0.6, 0.9),
                "HTR" => RGBf(0.9, 0.6, 0.2),
                "SFR" => RGBf(0.8, 0.3, 0.3)
            )

            type_color = get(type_colors, pj.type, RGBf(0.5, 0.5, 0.5))

            Box(fig[row, col],
                color = (type_color, 0.8),
                cornerradius = 3,
                strokewidth = 0,
                width = 40,
                height = 20,
                halign = :left,
                valign = :bottom,
                padding = (3, 3, 3, 3))

            Label(fig[row, col],
                  pj.type,
                  fontsize = 9,
                  font = :bold,
                  color = :white,
                  halign = :left,
                  valign = :bottom,
                  padding = (8, 8, 8, 8))

            idx += 1
        end
    end

    # Add shared legend at top
    legend_elements = [
        LineElement(color = lr_colors[0.05], linestyle = lr_styles[0.05], linewidth = 2),
        LineElement(color = lr_colors[0.10], linestyle = lr_styles[0.10], linewidth = 2),
        LineElement(color = lr_colors[0.15], linestyle = lr_styles[0.15], linewidth = 2),
        LineElement(color = :red, linestyle = :dashdot, linewidth = 2)
    ]

    legend_labels = [
        "LR 5%",
        "LR 10% (Base)",
        "LR 15%",
        "Vendor Baseline"
    ]

    Legend(fig[0, :],
           legend_elements,
           legend_labels,
           orientation = :horizontal,
           tellwidth = false,
           tellheight = true,
           framevisible = true,
           titlesize = 12,
           labelsize = 11)

    # Add reactor type legend
    type_legend_elements = [
        PolyElement(color = type_colors["BWR"]),
        PolyElement(color = type_colors["HTR"]),
        PolyElement(color = type_colors["PWR"]),
        PolyElement(color = type_colors["SFR"])
    ]

    type_legend_labels = ["BWR", "HTR", "PWR", "SFR"]

    Legend(fig[end+1, :],
           type_legend_elements,
           type_legend_labels,
           "Reactor Type Colors",
           orientation = :horizontal,
           tellwidth = false,
           tellheight = true,
           framevisible = true,
           titlesize = 11,
           labelsize = 10)

    return fig, crossing_data
end

"""
    learning_curves_pairs_appendix(pjs, wacc_range, electricity_price_mean, opt_scaling, construction_time_ranges)

Create multiple 2x1 (side-by-side) learning curve figures for appendix.
Instead of one large 2x5 grid, generates 5 separate figures with 2 reactors each.

# Arguments
- `pjs::Vector{project}`: Vector of reactor projects (will be processed in pairs)
- `wacc_range::Vector`: [min, max] WACC range
- `electricity_price_mean::Float64`: Mean electricity price for reference
- `opt_scaling::String`: Scaling method ("roulstone" or "rothwell")
- `construction_time_ranges::Dict`: Construction time ranges by scale

# Returns
- `figures::Vector{Figure}`: Vector of figures, each with 2 reactors side-by-side
- `crossing_data::DataFrame`: Crossing point data for all reactors
"""
function learning_curves_pairs_appendix(pjs::Vector, wacc_range::Vector,
                                       electricity_price_mean::Float64,
                                       opt_scaling::String,
                                       construction_time_ranges::Dict)

    @info "Creating learning curves as paired figures (2x1 layout per figure)"

    n_reactors = length(pjs)
    n_pairs = ceil(Int, n_reactors / 2)

    figures = Figure[]
    crossing_data = DataFrame(
        Reactor = String[],
        VendorBaseline = Float64[],
        CrossingN = Union{Int, Missing}[],
        CrossingLCOE = Union{Float64, Missing}[]
    )

    # Colors and styles (defined once for consistency)
    lr_colors = Dict(
        0.05 => RGBf(0.7, 0.7, 0.7),
        0.10 => RGBf(0.2, 0.4, 0.8),
        0.15 => RGBf(0.3, 0.7, 0.3)
    )

    lr_styles = Dict(
        0.05 => :dash,
        0.10 => :solid,
        0.15 => :dot
    )

    for pair_idx in 1:n_pairs
        # Get reactors for this pair
        start_idx = (pair_idx - 1) * 2 + 1
        end_idx = min(start_idx + 1, n_reactors)
        pair_pjs = pjs[start_idx:end_idx]

        @info "  Creating figure $(pair_idx)/$(n_pairs) with $(length(pair_pjs)) reactor(s)"

        # Create figure for this pair - wider format for 2 side-by-side plots
        fig = Figure(size=(1100, 500), backgroundcolor=:white)

        for (col, pj) in enumerate(pair_pjs)
            @info "    Processing $(pj.name)"

            # Calculate vendor baseline using OPTIMISTIC assumptions
            vendor_baseline = calculate_vendor_baseline_lcoe_simple(pj, 0.07, 3)

            # Run learning simulation
            n_sim = 10000
            N_max = 30
            learning_rates = [0.05, 0.10, 0.15]
            N_values = 1:N_max

            results = Dict{Float64, Vector{Float64}}()

            for LR in learning_rates
                lcoe_trajectory = Float64[]

                for N in N_values
                    rand_vars = gen_rand_vars(opt_scaling, n_sim, wacc_range,
                                             electricity_price_mean, pj;
                                             apply_learning=true,
                                             N_unit=N,
                                             LR=LR,
                                             kappa=1.0,
                                             apply_soak_discount=false,
                                             construction_time_range=construction_time_ranges[pj.scale],
                                             quiet=true)

                    disc_res = mc_run(n_sim, pj, rand_vars)
                    res = npv_lcoe(disc_res, decompose=false)

                    push!(lcoe_trajectory, median(res.lcoe))
                end

                results[LR] = lcoe_trajectory
            end

            # Create subplot - row 1 for plot, col for position
            ax = Axis(fig[1, col],
                     xlabel = "Cumulative Units (N)",
                     ylabel = col == 1 ? "LCOE [EUR₂₀₂₅/MWh]" : "",
                     title = "$(pj.name) ($(pj.type))",
                     titlesize = 18,
                     titlefont = :bold,
                     xlabelsize = 14,
                     ylabelsize = 14,
                     xticklabelsize = 12,
                     yticklabelsize = 12)

            # Plot curves
            for LR in learning_rates
                lines!(ax, collect(N_values), results[LR],
                       color = lr_colors[LR],
                       linestyle = lr_styles[LR],
                       linewidth = 2.5)
            end

            # Vendor baseline
            hlines!(ax, [vendor_baseline],
                    color = :red,
                    linestyle = :dashdot,
                    linewidth = 2.5)

            # Find crossing for base case (LR=10%)
            base_trajectory = results[0.10]
            crossing_idx = findfirst(x -> x <= vendor_baseline, base_trajectory)

            if !isnothing(crossing_idx)
                crossing_n = N_values[crossing_idx]
                crossing_lcoe = base_trajectory[crossing_idx]

                scatter!(ax, [crossing_n], [crossing_lcoe],
                        color = :red,
                        markersize = 12,
                        marker = :star5)

                text!(ax, crossing_n + 1, crossing_lcoe,
                      text = "N=$(crossing_n)",
                      align = (:left, :center),
                      fontsize = 12,
                      color = :red,
                      font = :bold)

                push!(crossing_data, (pj.name, vendor_baseline, crossing_n, crossing_lcoe))
            else
                push!(crossing_data, (pj.name, vendor_baseline, missing, missing))
            end
        end

        # Add shared legend at bottom
        legend_elements = [
            LineElement(color = lr_colors[0.05], linestyle = lr_styles[0.05], linewidth = 2.5),
            LineElement(color = lr_colors[0.10], linestyle = lr_styles[0.10], linewidth = 2.5),
            LineElement(color = lr_colors[0.15], linestyle = lr_styles[0.15], linewidth = 2.5),
            LineElement(color = :red, linestyle = :dashdot, linewidth = 2.5)
        ]

        legend_labels = [
            "LR 5% (Pessimistic)",
            "LR 10% (Base)",
            "LR 15% (Optimistic)",
            "Vendor Baseline"
        ]

        Legend(fig[2, :],
               legend_elements,
               legend_labels,
               orientation = :horizontal,
               tellwidth = false,
               tellheight = true,
               framevisible = true,
               labelsize = 14,
               patchsize = (25, 10))

        push!(figures, fig)
    end

    @info "Created $(length(figures)) paired learning curve figures"

    return figures, crossing_data
end
