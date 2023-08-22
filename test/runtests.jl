import GGDUtils: interp
import OMAS: h5i2imas
import Statistics: mean
using Test

function test_interp()
    ids = h5i2imas("$(@__DIR__)/../samples/edge_profiles.h5")
    path_to_data = ["edge_profiles", "ggd", 1, "electrons", "density", 1, "values"]
    path_to_opd = ["edge_profiles", "grid_ggd", 1, "space", 1, "objects_per_dimension"]
    get_electron_density = interp(ids, path_to_data, path_to_opd)
    chosen_index = 555
    nodes = ids.edge_profiles.grid_ggd[1].space[1].objects_per_dimension[3].object[chosen_index].nodes
    nodes_coords = [ids.edge_profiles.grid_ggd[1].space[1].objects_per_dimension[1].object[node].geometry for node in nodes]
    cell_center = mean(nodes_coords)
    grid_val = ids.edge_profiles.ggd[1].electrons.density[1].values[chosen_index]
    searched_val = get_electron_density(cell_center...)
    println("Grid Value: ", grid_val)
    println("Searched Value: ", searched_val)
    @assert(grid_val == searched_val)
    return true
end

@testset "GGDUtils" begin
    @test test_interp()
end