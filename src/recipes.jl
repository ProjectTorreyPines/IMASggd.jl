@recipe function f(space::OMAS.edge_profiles__grid_ggd___space)
    nodes = space.objects_per_dimension[1].object
    edges = space.objects_per_dimension[2].object
    legend --> false
    linewidth --> 0.2
    linecolor --> :black
    size --> [600, 900]
    xaxis --> "R [m]"
    yaxis --> "Z [m]"
    label_assigned = false
    for edge in edges
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
                [Tuple(nodes[edge.nodes[ii]].geometry) for ii in [1, 2]]
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
    xaxis --> "R [m]"
    yaxis --> "Z [m]"
    subset_edge_inds = []
    label_assigned = false
    if subset.element[1].object[1].dimension == 2
        for ele in subset.element
            for bnd_ind in cells[ele.object[1].index].boundary
                union!(subset_edge_inds, bnd_ind)
            end
        end
    elseif subset.element[1].object[1].dimension == 1
        subset_edge_inds = [ele.object[1].index for ele in subset.element]
    elseif subset.element[1].object[1].dimension == 0
        for ele in subset.element
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
                ([nodes[ele.object[1].index].geometry[1]], [nodes[ele.object[1].index].geometry[2]])
            end
        end
    end
    for edge in edges[subset_edge_inds]
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
                [Tuple(nodes[edge.nodes[ii]].geometry) for ii in [1, 2]]
            end
        end
    end
end