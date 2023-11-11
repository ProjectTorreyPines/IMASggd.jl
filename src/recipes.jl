using RecipesBase
using ColorSchemes: ColorSchemes

@recipe function f(space::OMAS.edge_profiles__grid_ggd___space)
    nodes = space.objects_per_dimension[1].object
    edges = space.objects_per_dimension[2].object
    legend --> false
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
    linewidth --> 0.2
    linecolor --> :black
    size --> [600, 900]
    xaxis --> "R / m"
    yaxis --> "Z / m"
    subset_edge_inds = []
    label_assigned = false
    if subset.element[1].object[1].dimension == 3
        for ele ∈ subset.element
            for bnd_ind ∈ cells[ele.object[1].index].boundary
                union!(subset_edge_inds, bnd_ind)
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
    size --> [600, 900]
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
