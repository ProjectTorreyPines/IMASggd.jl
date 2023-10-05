import NearestNeighbors: KDTree, knn
import StaticArrays: SVector
import Statistics: mean

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

function get_kdtree(
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
)
    subset_centers = get_subset_centers(space, subset)
    return KDTree([SVector{2}(sc) for sc ∈ subset_centers]; leafsize=10)
end

"""
    interp(
    prop_values::Vector{T},
    kdtree::KDTree;
    use_nearest_n=4,

) where {T <: Real}

Lowest level interpolation function. It takes a vector of property values and a KDTree
defined over a 2D space with the same number of nodes as the property values. It returns
a function that can be used to interpolate the property values at any point in the
space.
"""
function interp(
    prop_values::Vector{T},
    kdtree::KDTree;
    use_nearest_n::Int=4,
) where {T <: Real}
    function get_interp_val(x::Real, y::Real)
        nearest_indices, distances = knn(kdtree, Array([x, y]), use_nearest_n)
        values = [prop_values[ii] for ii ∈ nearest_indices]
        weights = 1 ./ distances
        if any(isinf.(weights))
            return values[distances.==0][1]
        end
        return sum(weights .* values) / sum(weights)
    end
    get_interp_val(xy::Tuple{Real, Real}) = get_interp_val(xy...)
    return get_interp_val
end

"""
    interp(
    prop_values::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space;
    use_nearest_n::Int=4,

) where {T <: Real}

If the whole space is provided instead of a kdtree, calculate the kdtree for whole
space. Again, here it is assumed that the property values are porvided for each node
of the space.
"""
function interp(
    prop_values::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space;
    use_nearest_n::Int=4,
) where {T <: Real}
    return interp(prop, get_kdtree(space); use_nearest_n)
end

"""
    interp(
    prop_values::Vector{Real},
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset;
    use_nearest_n::Int=4,

)

If a subset of the space is provided, calculate the kdtree for the subset. In this case
it is assumed that the property values are provided for each element of the subset.
"""
function interp(
    prop_values::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset;
    use_nearest_n::Int=4,
) where {T <: Real}
    return interp(prop_values, get_kdtree(space, subset); use_nearest_n)
end

"""
    interp(
    prop::edge_profiles__prop_on_subset,
    grid_ggd::OMAS.edge_profiles__grid_ggd,
    value_field::Symbol=:values;
    use_nearest_n::Int=4,

)

Example:
grid_ggd = dd.edge_profiles.grid_ggd[1]
get_electron_density = interp(dd.edge_profiles.ggd[1].electrons.density[1], grid_ggd)
get_e_field_par = interp(dd.edge_profiles.ggd[1].e_field[1], grid_ggd, :parallel)
"""
function interp(
    prop::edge_profiles__prop_on_subset,
    grid_ggd::OMAS.edge_profiles__grid_ggd,
    value_field::Symbol=:values;
    use_nearest_n::Int=4,
)
    subset = get_grid_subset_with_index(grid_ggd, prop.grid_subset_index)
    space = grid_ggd.space[subset.element[1].object[1].space]
    return interp(getfield(prop, value_field), space, subset; use_nearest_n)
end

"""
    interp(
    prop_arr::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    value_field::Symbol=:values;
    use_nearest_n::Int=4,

) where {T <: edge_profiles__prop_on_subset}

Example:
sol = get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 23)
get_electron_density = interp(dd.edge_profiles.ggd[1].electrons.density, space, sol)
"""
function interp(
    prop_arr::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    value_field::Symbol=:values;
    use_nearest_n::Int=4,
) where {T <: edge_profiles__prop_on_subset}
    prop = get_prop_with_grid_subset_index(prop_arr, subset.identifier.index)
    return interp(getfield(prop, value_field), space, subset; use_nearest_n)
end

"""
    interp(
    prop_arr::Vector{T},
    grid_ggd::OMAS.edge_profiles__grid_ggd,
    grid_subset_index::Int,
    value_field::Symbol=:values;
    use_nearest_n::Int=4,

) where {T <: edge_profiles__prop_on_subset}

Example:
get_n_e_sep = interp(dd.edge_profiles.ggd[1].electrons.density, grid_ggd, 16)
"""
function interp(
    prop_arr::Vector{T},
    grid_ggd::OMAS.edge_profiles__grid_ggd,
    grid_subset_index::Int,
    value_field::Symbol=:values;
    use_nearest_n::Int=4,
) where {T <: edge_profiles__prop_on_subset}
    prop = get_prop_with_grid_subset_index(prop_arr, grid_subset_index)
    subset = get_grid_subset_with_index(grid_ggd, grid_subset_index)
    space = grid_ggd.space[subset.element[1].object[1].space]
    return interp(getfield(prop, value_field), space, subset; use_nearest_n)
end
