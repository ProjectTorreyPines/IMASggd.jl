module GGDUtils

import NearestNeighbors: KDTree, knn
import StaticArrays: SVector
import Statistics: mean
import OMAS
import SOLPS2IMAS: get_subset_space
using RecipesBase

export interp
export get_kdtree
export project_prop_on_subset!


include("recipes.jl")


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


"""
    get_subset_centers(space::OMAS.edge_profiles__grid_ggd___space,
                            subset::OMAS.edge_profiles__grid_ggd___grid_subset)

Returns an array of tuples corresponding to (r,z) coordinates of the center of
cells or the center of edges in the subset space.
"""
function get_subset_centers(space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset)
    subset_space = get_subset_space(space, subset.element)
    grid_nodes = space.objects_per_dimension[1].object
    return [Tuple(mean([grid_nodes[node].geometry for node in obj.nodes])) for obj in subset_space]
end


"""
    project_prop_on_subset!(prop, from_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
                                 to_subset::OMAS.edge_profiles__grid_ggd___grid_subset;
                                 space::OMAS.edge_profiles__grid_ggd___space)

This function can be used to add another instance on a property vector representing the value
in a new subset that can be taken as a projection from an existing larger subset.
Arguments:
prop: A property like electrons.density that is a vector of objects with fields coefficients,
      grid_index, grid_subset_index, and values. The different instances in the vector correspond
      to different grid_subset for which the property is provided.
from_subset: grid_subset object which is already represented in the property instance. grid subset
             with index 5 is populated in electrons.density already if the values for all cells are
             present.
to_subset: grid_subset which is either a smaller part of from_subset (core, sol, idr, odr) but has
           same dimensions as from_subset
           OR
           is smaller in dimension that goes through the from_subset (core_boundary, separatix etc.)
space: (optional) space object in grid_ggd is required only when from_subset is higher dimensional
       than to_subset.
Returns:
NOTE: This function ends in ! which means it updates prop argument in place. But for the additional
utility, this function also returns a tuple
to_subset_centers, to_prop.values (when from_subset dimension is greater than to_subset dimension)
to_subset_ele_obj_inds, to_prop.values (when from_subset dimension is same as to_subset dimension)
Descriptions:
to_subset_centers: center of cells or center of edges of the to_subset where property values are
                   defined and stored
to_subset_ele_obj_inds: Indices of the elements of to_subset where property values are deined and
                        stored
to_prop.values: The projected values of the properties added to prop object in a new instance
"""
function project_prop_on_subset!(prop, from_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
                                 to_subset::OMAS.edge_profiles__grid_ggd___grid_subset;
                                 space::OMAS.edge_profiles__grid_ggd___space)
    if from_subset.element[1].object[1].dimension == to_subset.element[1].object[1].dimension
        return project_prop_on_subset!(prop, from_subset, to_subset)
    elseif from_subset.element[1].object[1].dimension > to_subset.element[1].object[1].dimension
        from_prop = nothing
        for p in prop
            if p.grid_subset_index == from_subset.identifier.index
                from_prop = p
                break
            end
        end
        if isnothing(from_prop)
            println("from_subset not represented in the property yet")
        end
        to_subset_centers = get_subset_centers(space, to_subset)
        resize!(prop, length(prop) + 1)
        to_prop = prop[end]
        to_prop.grid_index = from_prop.grid_index
        to_prop.grid_subset_index = to_subset.identifier.index
        resize!(to_prop.values, length(to_subset.element))
        prop_interp = interp(from_prop, space)
        to_prop.values = prop_interp(to_subset_centers)
        return to_subset_centers, to_prop.values
    else
        println("to_subset is higher dimensional than from_subset")
    end
end


function project_prop_on_subset!(prop, from_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
                                 to_subset::OMAS.edge_profiles__grid_ggd___grid_subset)
    from_prop = nothing
    for p in prop
        if p.grid_subset_index == from_subset.identifier.index
            from_prop = p
            break
        end
    end
    if isnothing(from_prop)
        println("from_subset not represented in the property yet")
    end
    if from_subset.element[1].object[1].dimension == to_subset.element[1].object[1].dimension
        resize!(prop, length(prop) + 1)
        to_prop = prop[end]
        to_prop.grid_index = from_prop.grid_index
        to_prop.grid_subset_index = to_subset.identifier.index
        from_subset_ele_obj_inds = [ele.object[1].index for ele in from_subset.element]
        to_subset_ele_obj_inds = [ele.object[1].index for ele in to_subset.element]
        if to_subset_ele_obj_inds ⊆ from_subset_ele_obj_inds
            from_ele_inds = []
            for to_ele_obj_ind in to_subset_ele_obj_inds
                for (from_ele_ind, from_ele_obj_ind) in enumerate(from_subset_ele_obj_inds)
                    if from_ele_obj_ind == to_ele_obj_ind
                        append!(from_ele_inds, from_ele_ind)
                    end
                end
            end
            filtered_values = [from_prop.values[from_ele_ind] for from_ele_ind in from_ele_inds]
            resize!(to_prop.values, length(filtered_values)) 
            to_prop.values = filtered_values
            return to_subset_ele_obj_inds, to_prop.values
        else
            error("to_subset does not lie entirely inside from_subset. Projection not possible.")
        end
    else
        error("Dimensions of from_subset and to_subset do not match. Provide keyword ",
              "argument space if you want to project to a smaller dimension as space ",
              "information is required for that. Use\n",
              "project_prop_on_subset!(prop, from_subset, to_subset; space=space)")
    end
end

end # module GGDUtils
