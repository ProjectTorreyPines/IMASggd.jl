module GGDUtils

import NearestNeighbors: KDTree, knn
import SOLPS2IMAS: path_to_obj
import StaticArrays: SVector
import Statistics: mean

export interp

function interp(ids, path_to_data, path_to_opd)
    prop = path_to_obj(ids, path_to_data)
    grid_nodes = path_to_obj(ids, path_to_opd)[1].object
    grid_faces = path_to_obj(ids, path_to_opd)[3].object
    grid_faces = [cell for cell in grid_faces if length(cell.nodes) == 4]
    grid_centers = [SVector{2}(mean([grid_nodes[node].geometry for node in cell.nodes])) for cell in grid_faces]
    kdtree = KDTree(grid_centers; leafsize=10)
    function get_interp_val(x, y)
        nearest_indices, distances = knn(kdtree, Array([x, y]), 4)
        v1, v2, v3, v4 = [prop[ii] for ii in nearest_indices]
        d1, d2, d3, d4 = distances
        return ((v1 * d2 * d3 * d4 + d1 * v2 * d3 * d4
                 + d1 * d2 * v3 * d4 + d1 * d2 * d3 * v4)
                / (d2 * d3 * d4 + d1 * d3 * d4
                   + d1 * d2 * d4 + d1 * d2 * d3))
    end
    return get_interp_val
end

end # module GGDUtils
