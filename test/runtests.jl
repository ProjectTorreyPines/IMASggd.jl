import GGDUtils: interp, get_kdtree, project_prop_on_subset!, get_grid_subset_with_index
import OMAS: h5i2imas
import Statistics: mean
import SOLPS2IMAS: solps2imas
using Test
using ArgParse: ArgParse

allowed_rtol = 1e-4

function parse_commandline()
    s = ArgParse.ArgParseSettings(; description="Run tests. Default is all tests.")

    ArgParse.add_arg_table!(s,
        ["--interp"],
        Dict(:help => "Test interp",
            :action => :store_true),
        ["--projection"],
        Dict(:help => "Test project_prop_on_subset!()",
            :action => :store_true),
        ["--in"],
        Dict(:help => "Test ∈",
            :action => :store_true),
        ["--interpeqt"],
        Dict(:help => "Test interpolation of equilibrium time slice",
            :action => :store_true),
    )
    args = ArgParse.parse_args(s)
    if !any(values(args)) # If no flags are set, run all tests
        for k ∈ keys(args)
            args[k] = true
        end
    end
    return args
end
args = parse_commandline()

if args["interp"]
    @testset "interp" begin
        # ids = h5i2imas("$(@__DIR__)/../samples/edge_profiles.h5")
        b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
        b2output = "$(@__DIR__)/../samples/b2time.nc"
        gsdesc = "$(@__DIR__)/../samples/gridspacedesc.yml"
        b2mn = "$(@__DIR__)/../samples/b2mn.dat"
        ids = solps2imas(b2gmtry, b2output, gsdesc, b2mn)
        n_e = ids.edge_profiles.ggd[1].electrons.density[1]
        grid_ggd = ids.edge_profiles.grid_ggd[1]
        space = grid_ggd.space[1]

        chosen_index = 555
        nodes = space.objects_per_dimension[3].object[chosen_index].nodes
        nodes_coords =
            [space.objects_per_dimension[1].object[node].geometry for node ∈ nodes]
        cell_center = mean(nodes_coords)
        grid_val = ids.edge_profiles.ggd[1].electrons.density[1].values[chosen_index]

        # test interp(prop, grid_ggd)
        get_n_e = interp(n_e, grid_ggd)
        searched_val = get_n_e(cell_center...)
        println("Electron density at: ", cell_center)
        println("Grid Value: ", grid_val)
        println("Searched Value: ", searched_val)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        # test interp(prop_arr, space, subset)
        subset = get_grid_subset_with_index(grid_ggd, 5)
        get_n_e = interp(ids.edge_profiles.ggd[1].electrons.density, space, subset)
        searched_val = get_n_e(cell_center...)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        # test interp(prop_arr, grid_ggd, grid_subset_index)
        get_n_e = interp(ids.edge_profiles.ggd[1].electrons.density, grid_ggd, 5)
        searched_val = get_n_e(cell_center...)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        # Use the kdtree to interpolate several quantities using
        # inverse distance weighing
        kdtree = get_kdtree(space)
        get_T_e =
            interp(ids.edge_profiles.ggd[1].electrons.temperature[1].values, kdtree)
        get_n_e = interp(ids.edge_profiles.ggd[1].electrons.density[1].values, kdtree)

        grid_val =
            ids.edge_profiles.ggd[1].electrons.temperature[1].values[chosen_index]
        searched_val = get_T_e(cell_center...)
        println("Electron temperature at: ", cell_center)
        println("Grid Value: ", grid_val)
        println("Searched Value: ", searched_val)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        grid_val = ids.edge_profiles.ggd[1].electrons.density[1].values[chosen_index]
        searched_val = get_n_e(cell_center...)
        println("Electron density at: ", cell_center)
        println("Grid Value: ", grid_val)
        println("Searched Value: ", searched_val)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        chosen_index = 553:557
        grid_val = ids.edge_profiles.ggd[1].electrons.density[1].values[chosen_index]
        chosen_nodes =
            [space.objects_per_dimension[3].object[ii].nodes for ii ∈ chosen_index]
        cell_centers = [
            Tuple(
                mean([
                    space.objects_per_dimension[1].object[node].geometry for
                    node ∈ nodes
                ]),
            ) for nodes ∈ chosen_nodes
        ]
        searched_val = get_n_e.(cell_centers)
        println("Electron density at: ", cell_center)
        println("Grid Value: ", grid_val)
        println("Searched Value: ", searched_val)
        @test mean(abs.((grid_val .- searched_val) ./ grid_val)) < allowed_rtol
    end
end

if args["projection"]
    @testset "project_prop_on_subset!" begin
        b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
        b2output = "$(@__DIR__)/../samples/b2time.nc"
        gsdesc = "$(@__DIR__)/../samples/gridspacedesc.yml"
        b2mn = "$(@__DIR__)/../samples/b2mn.dat"
        dd = solps2imas(b2gmtry, b2output, gsdesc, b2mn)
        space = dd.edge_profiles.grid_ggd[1].space[1]
        prop = dd.edge_profiles.ggd[1].electrons.density
        # All cells
        from_subset = get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 5)
        kdtree = get_kdtree(space, from_subset)
        # separatix
        to_subset = get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 16)
        separatix_centers, values_at_separatix =
            project_prop_on_subset!(prop, from_subset, to_subset, space)
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
        @test true
    end
end

if args["in"]
    @testset "test ∈" begin
        b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
        b2output = "$(@__DIR__)/../samples/b2time.nc"
        gsdesc = "$(@__DIR__)/../samples/gridspacedesc.yml"
        b2mn = "$(@__DIR__)/../samples/b2mn.dat"
        dd = solps2imas(b2gmtry, b2output, gsdesc, b2mn)
        grid_ggd = dd.edge_profiles.grid_ggd[1]
        space = grid_ggd.space[1]
        subset_corebnd = get_grid_subset_with_index(grid_ggd, 15)
        subset_sol = get_grid_subset_with_index(grid_ggd, 23)
        subset_odr = get_grid_subset_with_index(grid_ggd, 24)

        @test (6.0, 0.0) ∈ (subset_corebnd, space)
        @test (5.0, -2.5) ∉ (subset_corebnd, space)
        @test (6.0, 4.0) ∈ (subset_sol, space)
        @test (6.0, 3.0) ∉ (subset_sol, space)
        @test (5.1, -3.7) ∈ (subset_odr, space)
        @test (4.5, -3.7) ∉ (subset_odr, space)
    end
end
