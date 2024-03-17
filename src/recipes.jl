using RecipesBase
using ColorSchemes: ColorSchemes
import Statistics: norm, dot

@recipe function f(space::OMAS.edge_profiles__grid_ggd___space)
    nodes = space.objects_per_dimension[1].object
    edges = space.objects_per_dimension[2].object
    legend --> false
    subplot --> 1
    linewidth --> 0.2
    linecolor --> :black
    size --> [600, 900]
    xaxis --> "R / m"
    yaxis --> "Z / m"
    label_assigned = false
    for edge ∈ edges
        if 0 ∉ edge.nodes
            @series begin
                seriestype := :path
                if !label_assigned
                    if :label ∈ keys(plotattributes)
                        label := plotattributes[:label]
                    else
                        label := space.identifier.name
                    end
                    label_assigned = true
                else
                    label := ""
                end
                [Tuple(nodes[edge.nodes[ii]].geometry) for ii ∈ [1, 2]]
            end
        end
    end
end

@recipe function f(
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
)
    nodes = space.objects_per_dimension[1].object
    edges = space.objects_per_dimension[2].object
    cells = space.objects_per_dimension[3].object
    legend --> false
    subplot --> 1
    linewidth --> 0.2
    linecolor --> :black
    size --> [600, 900]
    xaxis --> "R / m"
    yaxis --> "Z / m"
    subset_edge_inds = []
    label_assigned = false
    if subset.element[1].object[1].dimension == 3
        for ele ∈ subset.element
            for bnd ∈ cells[ele.object[1].index].boundary
                union!(subset_edge_inds, bnd.index)
            end
        end
    elseif subset.element[1].object[1].dimension == 2
        subset_edge_inds = [ele.object[1].index for ele ∈ subset.element]
    elseif subset.element[1].object[1].dimension == 1
        for ele ∈ subset.element
            @series begin
                seriestype := :scatter
                marker --> (:circle, 5)
                if !label_assigned
                    if :label ∈ keys(plotattributes)
                        label := plotattributes[:label]
                    else
                        label := subset.identifier.name
                    end
                    label_assigned = true
                else
                    label := ""
                end
                (
                    [nodes[ele.object[1].index].geometry[1]],
                    [nodes[ele.object[1].index].geometry[2]],
                )
            end
        end
    end
    for edge ∈ edges[subset_edge_inds]
        if 0 ∉ edge.nodes
            @series begin
                seriestype := :path
                if !label_assigned
                    if :label ∈ keys(plotattributes)
                        label := plotattributes[:label]
                    else
                        label := subset.identifier.name
                    end
                    label_assigned = true
                else
                    label := ""
                end
                [Tuple(nodes[edge.nodes[ii]].geometry) for ii ∈ [1, 2]]
            end
        end
    end
end

@recipe function f(grid_ggd::OMAS.edge_profiles__grid_ggd, prop::OMAS.IDSvectorElement)
    subset = get_grid_subset_with_index(grid_ggd, prop.grid_subset_index)
    space = grid_ggd.space[subset.element[1].object[1].space]
    nodes = space.objects_per_dimension[1].object
    cells = space.objects_per_dimension[3].object
    legend --> false
    size --> [635, 900]
    xaxis --> "R / m"
    yaxis --> "Z / m"
    layout := @layout [a{0.95w} b]
    if subset.element[1].object[1].dimension == 3
        if :seriescolor in keys(plotattributes)
            color_scheme = plotattributes[:seriescolor]
        else
            color_scheme = :inferno
        end
        color_grad = getproperty(ColorSchemes, color_scheme)
        val_min = minimum(prop.values)
        val_max = maximum(prop.values)
        function get_color(prop_value)
            return color_grad[(log10(prop_value / val_min)/log10(val_max / val_min))]
        end

        if :colorbar_title in keys(plotattributes)
            prop_name = plotattributes[:colorbar_title]
        else
            prop_name = join(split(split(string(typeof(prop)), "___")[end], "__"), " ")
        end

        # Make a mock heatmap to get the colorbar in a side pane
        @series begin
            subplot := 2
            seriestype := :heatmap
            clims := (val_min, val_max)
            framestyle := :none
            color := color_scheme
            colorbar := true
            colorbar_scale := :log10
            colorbar_formatter := :scientific
            colorbar_title := prop_name
            lims := (-1, 0)
            rand(2, 2)
        end

        # Actual plot is plotted as different cell shapes filled with color values
        # based on property value
        for (ele, prop_value) ∈ zip(subset.element, prop.values)
            cell = cells[ele.object[1].index]
            @series begin
                subplot := 1
                seriestype := :shape
                linecolor := get_color(prop_value)
                fillcolor := get_color(prop_value)
                label := ""
                [Tuple(nodes[cell.nodes[ii]].geometry) for ii ∈ [1, 2, 4, 3]]
            end
        end
    end
end

@recipe function f(
    grid_ggd_arr::Vector{OMAS.edge_profiles__grid_ggd},
    prop::OMAS.IDSvectorElement,
)
    found = false
    for grid_ggd ∈ grid_ggd_arr
        if grid_ggd.identifier.index == prop.grid_index
            @series begin
                return grid_ggd, prop
            end
            found = true
        end
    end
    if !found
        error("Provided property belongs to a grid_index tha is not present in ",
            "provided grid_ggd array",
        )
    end
end

@recipe f(
    ifo::OMAS.interferometer,
) =
    for ch ∈ ifo.channel
        @series begin
            return ch
        end
    end

@recipe function f(
    ifo_ch::OMAS.interferometer__channel,
)
    if :plot_type ∈ keys(plotattributes)
        plot_type = plotattributes[:plot_type]
    else
        plot_type = :los
    end
    if plot_type == :los
        @series begin
            label --> ifo_ch.name
            ifo_ch.line_of_sight
        end
    elseif plot_type == :n_e || plot_type == :n_e_line
        @series begin
            label --> ifo_ch.name
            ifo_ch.n_e_line
        end
    elseif plot_type == :n_e_average || plot_type == :n_e_line_average
        @series begin
            label --> ifo_ch.name
            ifo_ch.n_e_line_average
        end
    else
        error("Invalid plot type. Choose between los (line of sight, R-Z plane), ",
            "n_e (integrated n_e along line of sight vs time), ",
            "and n_e_average (average n_e vs time)")
    end
end

@recipe function f(
    ifo_ch_los::OMAS.interferometer__channel___line_of_sight,
)
    subplot --> 1
    size --> [600, 900]
    xaxis --> "R / m"
    yaxis --> "Z / m"
    fp = ifo_ch_los.first_point
    sp = ifo_ch_los.second_point
    tp = ifo_ch_los.third_point
    if tp == OMAS.interferometer__channel___line_of_sight__third_point()
        tp = fp
    end

    if :mirror ∈ keys(plotattributes)
        mirror = plotattributes[:mirror]
    else
        mirror = true
    end
    if mirror
        # Calculate mirror points
        fsuv = Array([sp.r - fp.r, sp.z - fp.z])
        stuv = Array([tp.r - sp.r, tp.z - sp.z])
        fsuv = fsuv / norm(fsuv)
        stuv = stuv / norm(stuv)
        mirror_uv = fsuv + stuv
        if mirror_uv[1] == 0 && mirror_uv[2] == 0
            mirror_uv = [fsuv[2], -fsuv[1]]
        end
        mirror_uv = mirror_uv / norm(mirror_uv)
        if :mirror_length ∈ keys(plotattributes)
            mirror_length = plotattributes[:mirror_length]
        else
            mirror_length = 0.5
        end
        mirror_perp = ([-mirror_uv[2], mirror_uv[1]])
        mirror_perp = mirror_perp * sign(dot(mirror_perp, fsuv))
        if :mirror_thickness ∈ keys(plotattributes)
            mirror_thickness = plotattributes[:mirror_thickness]
        else
            mirror_thickness = 0.1
        end
        mirror_r1 = sp.r - mirror_uv[1] * mirror_length / 2
        mirror_r2 = sp.r + mirror_uv[1] * mirror_length / 2
        mirror_r3 = mirror_r2 + mirror_perp[1] * mirror_thickness
        mirror_r4 = mirror_r1 + mirror_perp[1] * mirror_thickness
        mirror_z1 = sp.z - mirror_uv[2] * mirror_length / 2
        mirror_z2 = sp.z + mirror_uv[2] * mirror_length / 2
        mirror_z3 = mirror_z2 + mirror_perp[2] * mirror_thickness
        mirror_z4 = mirror_z1 + mirror_perp[2] * mirror_thickness
    end

    # Draw line of sight
    @series begin
        seriestype := :path
        linewidth --> 2
        [fp.r, sp.r, tp.r], [fp.z, sp.z, tp.z]
    end

    if mirror
        # Draw mirror
        @series begin
            seriestype := :path
            label := ""
            linecolor := :gray
            linewidth := 2
            [mirror_r1, mirror_r2], [mirror_z1, mirror_z2]
        end
        @series begin
            seriestype := :shape
            label := ""
            linealpha := 0
            fillstyle := :/
            fillcolor := :black
            [mirror_r1, mirror_r2, mirror_r3, mirror_r4],
            [mirror_z1, mirror_z2, mirror_z3, mirror_z4]
        end
    end
end

@recipe function f(
    ifo_ch_n_e_line::OMAS.interferometer__channel___n_e_line,
)
    if :average ∈ keys(plotattributes)
        @series begin
            ifo_ch.n_e_line_average
        end
    else
        subplot --> 1
        xaxis --> "time / s"
        yaxis --> "Integrrated n_e / m^-2"
        @series begin
            seriestype := :path
            linewidth --> 2
            ifo_ch_n_e_line.time, ifo_ch_n_e_line.data
        end
    end
end

@recipe function f(
    ifo_ch_n_e_line_average::OMAS.interferometer__channel___n_e_line_average,
)
    subplot --> 1
    xaxis --> "time / s"
    yaxis --> "Average n_e / m^-3"
    @series begin
        seriestype := :path
        linewidth --> 2
        ifo_ch_n_e_line_average.time, ifo_ch_n_e_line_average.data
    end
end
