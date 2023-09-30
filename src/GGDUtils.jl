module GGDUtils

import NearestNeighbors: KDTree, knn
import StaticArrays: SVector
import Statistics: mean
using OMAS: OMAS
import SOLPS2IMAS: get_subset_space, get_grid_subset_with_index
using RecipesBase
using ColorSchemes: ColorSchemes

export interp
export get_kdtree
export project_prop_on_subset!
export get_subset_centers
export get_prop_with_grid_subset_index

include("recipes.jl")

function get_kdtree(space::OMAS.edge_profiles__grid_ggd___space)
    grid_nodes = space.objects_per_dimension[1].object
    grid_faces = space.objects_per_dimension[3].object
    grid_faces = [cell for cell ∈ grid_faces if length(cell.nodes) == 4]
    grid_centers = [
        SVector{2}(mean([grid_nodes[node].geometry for node ∈ cell.nodes])) for
        cell ∈ grid_faces
    ]
    return KDTree(grid_centers; leafsize=10)
end

function interp(prop, kdtree::KDTree)
    function get_interp_val(x::Number, y::Number)
        nearest_indices, distances = knn(kdtree, Array([x, y]), 4)
        v1, v2, v3, v4 = [prop.values[ii] for ii ∈ nearest_indices]
        d1, d2, d3, d4 = distances
        return (
            (
                v1 * d2 * d3 * d4 + d1 * v2 * d3 * d4
                + d1 * d2 * v3 * d4 + d1 * d2 * d3 * v4
            )
            /
            (d2 * d3 * d4 + d1 * d3 * d4
             + d1 * d2 * d4 + d1 * d2 * d3)
        )
    end
    function get_interp_val(x::Vector{Float64}, y::Vector{Float64})
        if length(x) != length(y)
            error("Length of the two axes are not equal")
        else
            return [get_interp_val(x[ii], y[ii]) for ii ∈ eachindex(x)]
        end
    end
    function get_interp_val(xy::Vector{Tuple{Float64, Float64}})
        return [get_interp_val(xy[ii]...) for ii ∈ eachindex(xy)]
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
    return [
        Tuple(mean([grid_nodes[node].geometry for node ∈ obj.nodes])) for
        obj ∈ subset_space
    ]
end

#! format: off
"""
    project_prop_on_subset!(prop,
    from_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    to_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    space::OMAS.edge_profiles__grid_ggd___space)

This function can be used to add another instance on a property vector representing the
value in a new subset that can be taken as a projection from an existing larger subset.
Arguments:
prop: A property like electrons.density that is a vector of objects with fields
      coefficients, grid_index, grid_subset_index, and values. The different instances
      in the vector correspond to different grid_subset for which the property is
      provided.
from_subset: grid_subset object which is already represented in the property instance.
             grid subset with index 5 is populated in electrons.density already if the
             values for all cells are present.
to_subset: grid_subset which is either a smaller part of from_subset (core, sol, idr,
           odr) but has same dimensions as from_subset
           OR
           is smaller in dimension that goes through the from_subset (core_boundary,
           separatix etc.)
space: (optional) space object in grid_ggd is required only when from_subset is
       higher dimensional than to_subset.
Returns:
NOTE: This function ends in ! which means it updates prop argument in place. But for
the additional utility, this function also returns a tuple
(to_subset_centers, to_prop.values) when from_subset dimension is greater than
                                    to_subset dimension
OR
(to_subset_ele_obj_inds, to_prop.values) when from_subset dimension is same as
                                         to_subset dimension)
Descriptions:
to_subset_centers: center of cells or center of edges of the to_subset where property
                   values are defined and stored
to_subset_ele_obj_inds: Indices of the elements of to_subset where property values are
                        defined and stored
to_prop.values: The projected values of the properties added to prop object in a new
                instance
"""
#! format: on
function project_prop_on_subset!(prop,
    from_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    to_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    space::OMAS.edge_profiles__grid_ggd___space)
    if from_subset.element[1].object[1].dimension ==
       to_subset.element[1].object[1].dimension
        return project_prop_on_subset!(prop, from_subset, to_subset)
    elseif from_subset.element[1].object[1].dimension >
           to_subset.element[1].object[1].dimension
        from_prop = nothing
        for p ∈ prop
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

function project_prop_on_subset!(prop,
    from_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    to_subset::OMAS.edge_profiles__grid_ggd___grid_subset,
)
    from_prop = get_prop_with_grid_subset_index(prop, from_subset.identifier.index)
    if isnothing(from_prop)
        println("from_subset not represented in the property yet")
    end
    if from_subset.element[1].object[1].dimension ==
       to_subset.element[1].object[1].dimension
        resize!(prop, length(prop) + 1)
        to_prop = prop[end]
        to_prop.grid_index = from_prop.grid_index
        to_prop.grid_subset_index = to_subset.identifier.index
        from_subset_ele_obj_inds = [ele.object[1].index for ele ∈ from_subset.element]
        to_subset_ele_obj_inds = [ele.object[1].index for ele ∈ to_subset.element]
        if to_subset_ele_obj_inds ⊆ from_subset_ele_obj_inds
            from_ele_inds = []
            for to_ele_obj_ind ∈ to_subset_ele_obj_inds
                for (from_ele_ind, from_ele_obj_ind) ∈
                    enumerate(from_subset_ele_obj_inds)
                    if from_ele_obj_ind == to_ele_obj_ind
                        append!(from_ele_inds, from_ele_ind)
                    end
                end
            end
            filtered_values =
                [from_prop.values[from_ele_ind] for from_ele_ind ∈ from_ele_inds]
            resize!(to_prop.values, length(filtered_values))
            to_prop.values = filtered_values
            return to_subset_ele_obj_inds, to_prop.values
        else
            error("to_subset does not lie entirely inside from_subset. Projection ",
                "not possible.",
            )
        end
    else
        error("Dimensions of from_subset and to_subset do not match. Provide keyword ",
            "argument space if you want to project to a smaller dimension as space ",
            "information is required for that. Use\n",
            "project_prop_on_subset!(prop, from_subset, to_subset; space=space)")
    end
end

function Base.:∈(
    point::Tuple{Float64, Float64},
    subset_of_space::Tuple{
        OMAS.edge_profiles__grid_ggd___grid_subset,
        OMAS.edge_profiles__grid_ggd___space,
    },
)
    r, z = point
    subset, space = subset_of_space
    dim = subset.element[1].object[1].dimension
    nodes = space.objects_per_dimension[1].object
    edges = space.objects_per_dimension[2].object
    if dim == 2
        subset_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
        subset_bnd.element = SOLPS2IMAS.get_subset_boundary(space, subset)
    elseif dim == 1
        subset_bnd = subset
    elseif dim == 0
        for ele ∈ subset.element
            node = nodes[ele.object[1].index]
            if node.geometry[1] == r && node.geometry[2] == z
                return true
            end
        end
        return false
    else
        error("Dimension ", dim, " is not supported yet.")
    end
    count = 0
    for ele ∈ subset_bnd.element
        edge = edges[ele.object[1].index]
        edge_z_at_r = line_between([nodes[node].geometry for node ∈ edge.nodes]...)(r)
        edge_ends =
            mapreduce(permutedims, vcat, [nodes[node].geometry for node ∈ edge.nodes])
        if maximum(edge_ends[:, 1]) >= r >= minimum(edge_ends[:, 1])
            if z < edge_z_at_r
                count += 1
            end
        end
    end
    if count % 2 == 1
        return true
    else
        return false
    end
end

function line_between(fp::Vector{Float64}, sp::Vector{Float64})
    slope = (sp[2] - fp[2]) / (sp[1] - fp[1])
    intercept = fp[2] - slope * fp[1]
    return (x::Float64) -> slope * x + intercept
end

function get_prop_with_grid_subset_index(prop, grid_subset_index::Int64)
    for p ∈ prop
        if p.grid_subset_index == grid_subset_index
            return p
        end
    end
    return nothing
end
end # module GGDUtils
