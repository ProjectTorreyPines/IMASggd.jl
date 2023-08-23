import GGDUtils: interp, get_kdtree
import OMAS: h5i2imas
import Statistics: mean
using Test

function test_interp()
    ids = h5i2imas("$(@__DIR__)/../samples/edge_profiles.h5")
    electron_density = ids.edge_profiles.ggd[1].electrons.density[1].values
    objects_per_dimension = ids.edge_profiles.grid_ggd[1].space[1].objects_per_dimension
    # path_to_data = ["edge_profiles", "ggd", 1, "electrons", "density", 1, "values"]
    # path_to_opd = ["edge_profiles", "grid_ggd", 1, "space", 1, "objects_per_dimension"]
    get_electron_density = interp(electron_density, objects_per_dimension)
    chosen_index = 555
    nodes = ids.edge_profiles.grid_ggd[1].space[1].objects_per_dimension[3].object[chosen_index].nodes
    nodes_coords = [ids.edge_profiles.grid_ggd[1].space[1].objects_per_dimension[1].object[node].geometry for node in nodes]
    cell_center = mean(nodes_coords)
    grid_val = ids.edge_profiles.ggd[1].electrons.density[1].values[chosen_index]
    searched_val = get_electron_density(cell_center...)
    println("Electron density at: ", cell_center)
    println("Grid Value: ", grid_val)
    println("Searched Value: ", searched_val)
    @assert(grid_val == searched_val)

    # Use the same kdtree to interpolate other quantities
    kdtree = get_kdtree(objects_per_dimension)
    get_electron_temperature = interp(ids.edge_profiles.ggd[1].electrons.temperature[1].values, kdtree)
    get_ion_density = interp(ids.edge_profiles.ggd[1].ion[1].density[1].values, kdtree)

    grid_val = ids.edge_profiles.ggd[1].electrons.temperature[1].values[chosen_index]
    searched_val = get_electron_temperature(cell_center...)
    println("Electron temperature at: ", cell_center)
    println("Grid Value: ", grid_val)
    println("Searched Value: ", searched_val)
    @assert(grid_val == searched_val)

    grid_val = ids.edge_profiles.ggd[1].ion[1].density[1].values[chosen_index]
    searched_val = get_ion_density(cell_center...)
    println("First Ion density at: ", cell_center)
    println("Grid Value: ", grid_val)
    println("Searched Value: ", searched_val)
    @assert(grid_val == searched_val)
    
    return true
end

@testset "GGDUtils" begin
    @test test_interp()
end