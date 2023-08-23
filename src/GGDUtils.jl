module GGDUtils

import NearestNeighbors: KDTree, knn
import StaticArrays: SVector
import Statistics: mean

export interp
export get_kdtree


function get_kdtree(objects_per_dimension)
    grid_nodes = objects_per_dimension[1].object
    grid_faces = objects_per_dimension[3].object
    grid_faces = [cell for cell in grid_faces if length(cell.nodes) == 4]
    grid_centers = [SVector{2}(mean([grid_nodes[node].geometry for node in cell.nodes])) for cell in grid_faces]
    return KDTree(grid_centers; leafsize=10)
end


function interp(prop, kdtree::KDTree)
    function get_interp_val(x, y)
        nearest_indices, distances = knn(kdtree, Array([x, y]), 4)
        v1, v2, v3, v4 = [prop[ii] for ii in nearest_indices]
        d1, d2, d3, d4 = distances
        return ((v1 * d2 * d3 * d4 + d1 * v2 * d3 * d4
                 + d1 * d2 * v3 * d4 + d1 * d2 * d3 * v4)
                /
                (d2 * d3 * d4 + d1 * d3 * d4
                 + d1 * d2 * d4 + d1 * d2 * d3))
    end
    return get_interp_val
end


function interp(prop, objects_per_dimension)
    return interp(prop, get_kdtree(objects_per_dimension))
end

end # module GGDUtils
