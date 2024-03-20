"""
    add_subset_element!(
    subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    sn::Int,
    dim::Int,
    index::Int,
    in_subset=(x...) -> true;
    kwargs...,

)

Adds a new element to gird_subset with properties space number (sn), dimension (dim),
and index (index). The element is added only if the function in_subset returns true.
"""
function add_subset_element!(
    subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    sn::Int,
    dim::Int,
    index::Int,
    in_subset=(x...) -> true;
    kwargs...,
)
    if in_subset(; kwargs...)
        dd_ind = length(subset.element) + 1
        resize!(subset.element, dd_ind)
        resize!(subset.element[dd_ind].object, 1)
        subset.element[dd_ind].object[1].space = sn
        subset.element[dd_ind].object[1].dimension = dim
        subset.element[dd_ind].object[1].index = index
    end
end

"""
    add_subset_element!(
    subset,
    sn,
    dim,
    index::Vector{Int},
    in_subset=(x...) -> true;
    kwargs...,

)

Overloaded to work differently (faster) with list of indices to be added.
"""
function add_subset_element!(
    subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    sn::Int,
    dim::Int,
    index::Vector{Int},
    in_subset=(x...) -> true;
    kwargs...,
)
    if in_subset(; kwargs...)
        dd_start_ind = length(subset.element) + 1
        resize!(subset.element, length(subset.element) + length(index))
        dd_stop_ind = length(subset.element)
        for (ii, dd_ind) ∈ enumerate(dd_start_ind:dd_stop_ind)
            resize!(subset.element[dd_ind].object, 1)
            subset.element[dd_ind].object[1].space = sn
            subset.element[dd_ind].object[1].dimension = dim
            subset.element[dd_ind].object[1].index = index[ii]
        end
    end
end

"""
    get_subset_space(space::IMASDD.edge_profiles__grid_ggd___space,
    elements::AbstractVector{<:IMASDD.edge_profiles__grid_ggd___grid_subset___element})

Returns an array of space object indices corresponding to the correct
objects_per_dimension (nodes, edges or cells) for the subset elements.
"""
function get_subset_space(space::IMASDD.edge_profiles__grid_ggd___space,
    elements::AbstractVector{<:IMASDD.edge_profiles__grid_ggd___grid_subset___element})
    nD = elements[1].object[1].dimension
    nD_objects = space.objects_per_dimension[nD].object
    return [nD_objects[ele.object[1].index] for ele ∈ elements]
end

"""
    get_grid_subset(
    grid_ggd::IMASDD.edge_profiles__grid_ggd,
    grid_subset_index::Int,

)

Returns the grid_subset in a grid_ggd with the matching grid_subset_index
"""
function get_grid_subset(
    grid_ggd::IMASDD.edge_profiles__grid_ggd,
    grid_subset_index::Int,
)
    for subset ∈ grid_ggd.grid_subset
        if subset.identifier.index == grid_subset_index
            return subset
        end
    end

    return error("Subset ", grid_subset_index, " not found.")
end

"""
    get_grid_subset(
    grid_ggd::IMASDD.edge_profiles__grid_ggd,
    grid_subset_name::String,

)

Returns the grid_subset in a grid_ggd with the matching grid_subset_name
"""
function get_grid_subset(
    grid_ggd::IMASDD.edge_profiles__grid_ggd,
    grid_subset_name::String,
)
    for subset ∈ grid_ggd.grid_subset
        if subset.identifier.name == grid_subset_name
            return subset
        end
    end

    return error("Subset ", grid_subset_name, " not found.")
end

"""
    get_subset_boundary_inds(
    space::IMASDD.edge_profiles__grid_ggd___space,
    subset::IMASDD.edge_profiles__grid_ggd___grid_subset,

)

Returns an array of space object indices corresponding to the boundary of the subset.
That means, it returns indices of nodes that are at the end of open edge subset or
it returns the indices of edges that are the the boundary of a cell subset.
Returns an empty array if the subset is 1D (nodes).
"""
function get_subset_boundary_inds(
    space::IMASDD.edge_profiles__grid_ggd___space,
    subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
)
    nD = subset.element[1].object[1].dimension
    if nD > 1  # Only 2D (edges) and 3D (cells) subsets have boundaries
        nD_objects = space.objects_per_dimension[nD].object
        elements = [nD_objects[ele.object[1].index] for ele ∈ subset.element]
        boundary_inds = Int[]
        for ele ∈ elements
            symdiff!(boundary_inds, [bnd.index for bnd ∈ ele.boundary])
        end
        return boundary_inds
    end
    return [] # 1D (nodes) subsets have no boundary
end

"""
    get_subset_boundary(
    space::IMASDD.edge_profiles__grid_ggd___space,
    subset::IMASDD.edge_profiles__grid_ggd___grid_subset,

)

Returns an array of elements of grid_subset generated from the boundary of the subset
provided. The dimension of these elments is reduced by 1.
"""
function get_subset_boundary(
    space::IMASDD.edge_profiles__grid_ggd___space,
    subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
)
    ret_subset = IMASDD.edge_profiles__grid_ggd___grid_subset()
    boundary_inds = get_subset_boundary_inds(space, subset)
    bnd_dim = subset.element[1].object[1].dimension - 1
    space_number = subset.element[1].object[1].space
    add_subset_element!(ret_subset, space_number, bnd_dim, boundary_inds)
    return ret_subset.element
end

"""
    subset_do(set_operator,
    itrs::Vararg{AbstractVector{<:IMASDD.edge_profiles__grid_ggd___grid_subset___element}};
    space::IMASDD.edge_profiles__grid_ggd___space=IMASDD.edge_profiles__grid_ggd___space(),
    use_nodes=false)

Function to perform any set operation (intersect, union, setdiff etc.) on
subset.element to generate a list of elements to go to subset object. If use_nodes is
true, the set operation will be applied on the set of nodes from subset.element, space
argument is required for this.
Note: that the arguments are subset.element (not the subset itself). Similarly, the
return object is a list of IMASDD.edge_profiles__grid_ggd___grid_subset___element.
"""
function subset_do(set_operator,
    itrs::Vararg{
        AbstractVector{<:IMASDD.edge_profiles__grid_ggd___grid_subset___element},
    };
    space::IMASDD.edge_profiles__grid_ggd___space=IMASDD.edge_profiles__grid_ggd___space(),
    use_nodes=false)
    if use_nodes
        ele_inds = set_operator(
            [
                union([
                        obj.nodes for obj ∈ get_subset_space(space, set_elements)
                    ]...
                ) for set_elements ∈ itrs
            ]...,
        )
        dim = 1
    else
        ele_inds = set_operator(
            [[ele.object[1].index for ele ∈ set_elements] for set_elements ∈ itrs]...,
        )
        dim = itrs[1][1].object[1].dimension
    end
    ret_subset = IMASDD.edge_profiles__grid_ggd___grid_subset()
    space_number = itrs[1][1].object[1].space
    add_subset_element!(ret_subset, space_number, dim, ele_inds)
    return ret_subset.element
end

"""
    get_subset_centers(space::IMASDD.edge_profiles__grid_ggd___space,
                            subset::IMASDD.edge_profiles__grid_ggd___grid_subset)

Returns an array of tuples corresponding to (r,z) coordinates of the center of
cells or the center of edges in the subset space.
"""
function get_subset_centers(space::IMASDD.edge_profiles__grid_ggd___space,
    subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
)
    subset_space = get_subset_space(space, subset.element)
    grid_nodes = space.objects_per_dimension[1].object
    return [
        Tuple(mean(SVector{2}(grid_nodes[node].geometry) for node ∈ obj.nodes)) for
        obj ∈ subset_space
    ]
end

#! format: off
"""
    project_prop_on_subset!(prop_arr::AbstractVector{T},
    from_subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    to_subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    space::IMASDD.edge_profiles__grid_ggd___space,
)

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
(to_subset_centers, to_prop_values) when from_subset dimension is greater than
                                    to_subset dimension
OR
(to_subset_ele_obj_inds, to_prop_values) when from_subset dimension is same as
                                         to_subset dimension)
Descriptions:
to_subset_centers: center of cells or center of edges of the to_subset where property
                   values are defined and stored
to_subset_ele_obj_inds: Indices of the elements of to_subset where property values are
                        defined and stored
to_prop_values: The projected values of the properties added to prop object in a new
                instance
"""
#! format: on
function project_prop_on_subset!(prop_arr::AbstractVector{T},
    from_subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    to_subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    space::IMASDD.edge_profiles__grid_ggd___space,
    value_field::Symbol=:values;
    interp_method::Symbol=:thin_plate_spline,
    interp_kwargs=Dict(),
) where {T <: edge_profiles__prop_on_subset}
    if from_subset.element[1].object[1].dimension ==
       to_subset.element[1].object[1].dimension
        return project_prop_on_subset!(prop, from_subset, to_subset)
    elseif from_subset.element[1].object[1].dimension >
           to_subset.element[1].object[1].dimension
        from_prop =
            get_prop_with_grid_subset_index(prop_arr, from_subset.identifier.index)
        if isnothing(from_prop)
            error("from_subset not represented in the property yet")
        end
        to_subset_centers = get_subset_centers(space, to_subset)
        resize!(prop_arr, length(prop_arr) + 1)
        to_prop = prop_arr[end]
        to_prop.grid_index = from_prop.grid_index
        to_prop.grid_subset_index = to_subset.identifier.index
        to_prop_values = getfield(to_prop, value_field)
        from_prop_values = getfield(from_prop, value_field)
        resize!(to_prop_values, length(to_subset.element))
        if interp_method == :thin_plate_spline
            prop_interp = interp(prop_arr, space, from_subset)
        elseif interp_method == :KDTree
            prop_interp = interp(
                from_prop_values,
                get_kdtree(space, from_subset; interp_kwargs...),
            )
        else
            error("Supported interpolation methods are :thin_plate_spline and :KDTree")
        end
        to_prop_values = prop_interp.(to_subset_centers)
        return to_subset_centers, to_prop_values
    else
        error("to_subset is higher dimensional than from_subset")
    end
end

function project_prop_on_subset!(prop_arr::AbstractVector{T},
    from_subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    to_subset::IMASDD.edge_profiles__grid_ggd___grid_subset,
    value_field::Symbol=:values,
) where {T <: edge_profiles__prop_on_subset}
    from_prop = get_prop_with_grid_subset_index(prop_arr, from_subset.identifier.index)
    if isnothing(from_prop)
        error("from_subset not represented in the property yet")
    end
    if from_subset.element[1].object[1].dimension ==
       to_subset.element[1].object[1].dimension
        resize!(prop_arr, length(prop_arr) + 1)
        to_prop = prop_arr[end]
        to_prop.grid_index = from_prop.grid_index
        to_prop.grid_subset_index = to_subset.identifier.index
        to_prop_values = getfield(to_prop, value_field)
        from_prop_values = getfield(from_prop, value_field)
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
                [from_prop_values[from_ele_ind] for from_ele_ind ∈ from_ele_inds]
            resize!(to_prop_values, length(filtered_values))
            to_prop_values = filtered_values
            return to_subset_ele_obj_inds, to_prop_values
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

"""
    Base.:∈(
    point::Tuple{Real, Real},
    subset_of_space::Tuple{
        IMASDD.edge_profiles__grid_ggd___grid_subset,
        IMASDD.edge_profiles__grid_ggd___space,
    },

)

Overloading ∈ operator to check if a point is inside a subset of space.

If the subset is 1-dimensional, all points are searched. If the subset is 2-dimensional,
it is checked if the point is within the enclosed area. It is assumed that a
2-dimensional subset used in such a context will form a closed area. If the subset is
3-dimensional, its boundary is calculated on the fly. If used multiple times, it is
recommended to calculate the boundary once and store it in a variable.
"""
function Base.:∈(
    point::Tuple{Real, Real},
    subset_of_space::Tuple{
        IMASDD.edge_profiles__grid_ggd___grid_subset,
        IMASDD.edge_profiles__grid_ggd___space,
    },
)
    r, z = point
    subset, space = subset_of_space
    dim = subset.element[1].object[1].dimension
    nodes = space.objects_per_dimension[1].object
    edges = space.objects_per_dimension[2].object
    if dim == 3
        subset_bnd = IMASDD.edge_profiles__grid_ggd___grid_subset()
        subset_bnd.element = get_subset_boundary(space, subset)
    elseif dim == 2
        subset_bnd = subset
    elseif dim == 1
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
    # Count number of times an upward going ray from (r,z) intersects the boundary
    count = 0
    for ele ∈ subset_bnd.element
        edge = edges[ele.object[1].index]
        r_max = maximum(nodes[node].geometry[1] for node ∈ edge.nodes)
        r_min = minimum(nodes[node].geometry[1] for node ∈ edge.nodes)
        if r_min <= r < r_max
            z_max = maximum(nodes[node].geometry[2] for node ∈ edge.nodes)
            if z < z_max
                count += 1
            end
        end
    end
    # If it is even, the point is outside the boundary
    return count % 2 == 1
end

function get_prop_with_grid_subset_index(
    prop::AbstractVector{T},
    grid_subset_index::Int,
) where {T <: edge_profiles__prop_on_subset}
    for p ∈ prop
        if p.grid_subset_index == grid_subset_index
            return p
        end
    end
    return nothing
end
