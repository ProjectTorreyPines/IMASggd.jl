module GGDUtils

import NearestNeighbors: KDTree, knn
import StaticArrays: SVector
import Statistics: mean
import OMAS
import SOLPS2IMAS: get_subset_space

export interp
export get_kdtree
export project_prop_on_subset!


function get_kdtree(space::OMAS.edge_profiles__grid_ggd___space)
    grid_nodes = space.objects_per_dimension[1].object
    grid_faces = space.objects_per_dimension[3].object
    grid_faces = [cell for cell in grid_faces if length(cell.nodes) == 4]
    grid_centers = [SVector{2}(mean([grid_nodes[node].geometry for node in cell.nodes])) for cell in grid_faces]
    return KDTree(grid_centers; leafsize=10)
end


function interp(prop, kdtree::KDTree)
    function get_interp_val(x::Number, y::Number)
        nearest_indices, distances = knn(kdtree, Array([x, y]), 4)
        v1, v2, v3, v4 = [prop.values[ii] for ii in nearest_indices]
        d1, d2, d3, d4 = distances
        return ((v1 * d2 * d3 * d4 + d1 * v2 * d3 * d4
                 + d1 * d2 * v3 * d4 + d1 * d2 * d3 * v4)
                /
                (d2 * d3 * d4 + d1 * d3 * d4
                 + d1 * d2 * d4 + d1 * d2 * d3))
    end
    function get_interp_val(x::Vector{Float64}, y::Vector{Float64})
        if length(x) != length(y)
            error("Length of the two axes are not equal")
        else
            return [get_interp_val(x[ii], y[ii]) for ii in eachindex(x)]
        end
    end
    function get_interp_val(xy::Vector{Tuple{Float64, Float64}})
        return [get_interp_val(xy[ii]...) for ii in eachindex(xy)]
    end
    return get_interp_val
end


function interp(prop, space::OMAS.edge_profiles__grid_ggd___space)
    return interp(prop, get_kdtree(space))
end


function project_prop_on_subset!(prop, from_subset, to_subset;
                                 space::OMAS.edge_profiles__grid_ggd___space)
    from_prop = nothing
    for p in prop
        if p.grid_subset_index == from_subset.identifier.index
            from_prop = p
            break
        end
    end
    if isnothing(from_prop)
        println("from_subset not represented in the property yet")
        return false
    end
    
    if from_subset.element[1].object[1].dimension == to_subset.element[1].object[1].dimension
        resize!(prop, length(prop) + 1)
        to_prop = prop[end]
        to_prop.grid_index = from_prop.grid_index
        to_prop.grid_subset_index = to_subset.identifier.index
        from_subset_element_inds = [ele.index for ele in from_subset.element]
        filtered_values = [from_prop[ele.index] for ele in to_subset.element if ele.index ∈ from_subset_element_inds]
        resize!(to_prop.values, length(filtered_values)) 
        to_prop.values = filtered_values
    elseif from_subset.element[1].object[1].dimension > to_subset.element[1].object[1].dimension
        resize!(prop, length(prop) + 1)
        to_prop = prop[end]
        to_prop.grid_index = from_prop.grid_index
        to_prop.grid_subset_index = to_subset.identifier.index
        resize!(to_prop.values, length(to_subset.element))
        prop_interp = interp(from_prop, space)
        grid_nodes = space.objects_per_dimension[1].object
        to_subset_space = get_subset_space(space, to_subset.element)
        to_prop.values = prop_interp([Tuple(mean([grid_nodes[node].geometry for node in obj.nodes])) for obj in to_subset_space])
    end
end

end # module GGDUtils
