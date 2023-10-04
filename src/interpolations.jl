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
