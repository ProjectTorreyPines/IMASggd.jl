import GGDUtils: interp, get_kdtree, project_prop_on_subset!
import OMAS: h5i2imas
import Statistics: mean
import SOLPS2IMAS: solps2imas, get_grid_subset_with_index
using Test

function test_interp()
    ids = h5i2imas("$(@__DIR__)/../samples/edge_profiles.h5")
    electron_density = ids.edge_profiles.ggd[1].electrons.density[1]
    space = ids.edge_profiles.grid_ggd[1].space[1]
    get_electron_density = interp(electron_density, space)
    chosen_index = 555
    nodes = space.objects_per_dimension[3].object[chosen_index].nodes
    nodes_coords =
        [space.objects_per_dimension[1].object[node].geometry for node ∈ nodes]
    cell_center = mean(nodes_coords)
    grid_val = ids.edge_profiles.ggd[1].electrons.density[1].values[chosen_index]
    searched_val = get_electron_density(cell_center...)
    println("Electron density at: ", cell_center)
    println("Grid Value: ", grid_val)
    println("Searched Value: ", searched_val)
    @assert(grid_val == searched_val)

    # Use the same kdtree to interpolate other quantities
    kdtree = get_kdtree(space)
    get_electron_temperature =
        interp(ids.edge_profiles.ggd[1].electrons.temperature[1], kdtree)
    get_ion_density = interp(ids.edge_profiles.ggd[1].ion[1].density[1], kdtree)

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

    chosen_index = 553:557
    grid_val = ids.edge_profiles.ggd[1].ion[1].density[1].values[chosen_index]
    chosen_nodes =
        [space.objects_per_dimension[3].object[ii].nodes for ii ∈ chosen_index]
    cell_centers = [
        Tuple(
            mean([
                space.objects_per_dimension[1].object[node].geometry for node ∈ nodes
            ])
        ) for nodes ∈ chosen_nodes
    ]
    searched_val = get_ion_density(cell_centers)
    return true
end


function test_project_prop_on_subset()
    b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
    b2output = "$(@__DIR__)/../samples/b2time_red.nc"
    gsdesc = "$(@__DIR__)/../samples/gridspacedesc.yml"
    b2mn = "$(@__DIR__)/../samples/b2mn.dat"
    dd = solps2imas(b2gmtry, b2output, gsdesc, b2mn)
    space = dd.edge_profiles.grid_ggd[1].space[1]
    prop = dd.edge_profiles.ggd[1].electrons.density
    # All cells
    from_subset = get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 5)
    # separatix
    to_subset = get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 16)
    # Cheating: need to fix solps2imas, prop did not get correct grid_subset_index
    prop[5].grid_subset_index = 5
    separatix_centers, values_at_separatix =
        project_prop_on_subset!(prop, from_subset, to_subset; space)
    # println("Projected to separatix:")
    # for ii ∈ eachindex(separatix_centers)
    #     println(separatix_centers[ii], ": ", values_at_separatix[ii])
    # end

    subset_core = get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 22)
    core_element_inds, values_at_core =
        project_prop_on_subset!(prop, from_subset, subset_core)
    # println("Project to core:")
    # for ii ∈ eachindex(core_element_inds)
    #     println(
    #         "Element index: ",
    #         core_element_inds[ii],
    #         " has value : ",
    #         values_at_core[ii]
    #     )
    # end
    return true
end

@testset "GGDUtils" begin
    @test test_interp()
    @test test_project_prop_on_subset()
end
