import IMASggd:
    interp, get_kdtree, project_prop_on_subset!, get_grid_subset, get_grid_ggd,
    get_subset_boundary, subset_do, deepcopy_subset, get_TPS_mats, get_space, IMASdd,
    mean
using Test

allowed_rtol = 1e-4

print("json2imas() time: ")
@time ids = IMASdd.json2imas(
    "$(@__DIR__)/../samples/time_dep_edge_profiles_last_step_only.json",
)

if isempty(ARGS) || "interp" in ARGS
    @testset "interp" begin
        b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
        b2output = "$(@__DIR__)/../samples/b2time.nc"
        gsdesc = "$(@__DIR__)/../samples/gridspacedesc.yml"
        b2mn = "$(@__DIR__)/../samples/b2mn.dat"
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
        print("interp(prop, grid_ggd) time: ")
        @time get_n_e = interp(n_e, grid_ggd)
        searched_val = get_n_e(cell_center...)
        println("Electron density at: ", cell_center)
        println("Grid Value: ", grid_val)
        println("Searched Value: ", searched_val)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        # test interp(prop_arr, space, subset)
        subset = get_grid_subset(grid_ggd, -5)
        print("interp(prop_arr, space, subset) time: ")
        @time get_n_e =
            interp(ids.edge_profiles.ggd[1].electrons.density, space, subset)
        searched_val = get_n_e(cell_center...)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        # test interp(prop_arr, grid_ggd, grid_subset_index)
        print("interp(prop_arr, grid_ggd, grid_subset_index) time: ")
        @time get_n_e = interp(ids.edge_profiles.ggd[1].electrons.density, grid_ggd, -5)
        searched_val = get_n_e(cell_center...)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        # Use the TPS_mats to interpolate several quantities using
        print("get_TPS_mats(space, subset) time: ")
        @time TPS_mats = get_TPS_mats(space, subset)
        print(
            "interp(prop_arr(for n_e), TPS_mats, grid_subset_index, value_field) time: ",
        )
        @time get_n_e =
            interp(ids.edge_profiles.ggd[1].electrons.density, TPS_mats, -5, :values)
        searched_val = get_n_e(cell_center...)
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol
        print(
            "interp(prop_arr(for T_e), TPS_mats, grid_subset_index, value_field) time: ",
        )
        @time get_T_e =
            interp(
                ids.edge_profiles.ggd[1].electrons.temperature,
                TPS_mats,
                -5,
                :values,
            )
        searched_val = get_T_e(cell_center...)
        grid_val =
            ids.edge_profiles.ggd[1].electrons.temperature[1].values[chosen_index]
        @test abs.((grid_val .- searched_val) ./ grid_val) < allowed_rtol

        # Use the kdtree to interpolate several quantities using
        # inverse distance weighing
        kdtree = get_kdtree(space)
        print("interp(prop_arr, kdtree) time: ")
        @time get_T_e =
            interp(ids.edge_profiles.ggd[1].electrons.temperature[1].values, kdtree)
        @time get_n_e =
            interp(ids.edge_profiles.ggd[1].electrons.density[1].values, kdtree)

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

if isempty(ARGS) || "projection" in ARGS
    @testset "Test project_prop_on_subset!" begin
        prop = ids.edge_profiles.ggd[1].electrons.density
        # All cells
        from_subset = get_grid_subset(ids.edge_profiles.grid_ggd[1], -5)
        # separatix
        to_subset = get_grid_subset(ids.edge_profiles.grid_ggd[1], 16)
        print("project_prop_on_subset!(prop, from_subset, to_subset) time: ")
        @time separatix_centers, values_at_separatix =
            project_prop_on_subset!(prop, from_subset, to_subset)
        # println("Projected to separatix:")
        # for ii ∈ eachindex(separatix_centers)
        #     println(separatix_centers[ii], ": ", values_at_separatix[ii])
        # end

        subset_core =
            get_grid_subset(ids.edge_profiles.grid_ggd[1], 22)
        print("project_prop_on_subset!(prop, from_subset, subset_core) time: ")
        @time core_element_inds, values_at_core =
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

        idstd = IMASdd.json2imas(
            "$(@__DIR__)/../samples/time_dep_edge_profiles_with_interferometer.json",
        )
        # All cells
        from_subset = get_grid_subset(idstd.edge_profiles.grid_ggd[1], -5)
        # separatix
        to_subset = get_grid_subset(idstd.edge_profiles.grid_ggd[1], 16)
        print("project_prop_on_subset!(ggds, prop_path, from_subset, to_subset) time: ")
        @time projection_return =
            project_prop_on_subset!(
                idstd.edge_profiles.ggd,
                "electrons.density",
                from_subset,
                to_subset,
            )
        @test length(projection_return) == length(idstd.edge_profiles.ggd)
    end
end

if isempty(ARGS) || "subset_tools" in ARGS
    @testset "Test subset tools" begin
        grid_ggd = ids.edge_profiles.grid_ggd[1]
        space = grid_ggd.space[1]

        print("get_grid_ggd(space) time: ")
        @time grid_ggd_copy = get_grid_ggd(space)
        @test grid_ggd == grid_ggd_copy

        print("get_grid_subset(grid_ggd, 22) time: ")
        @time subset_core = get_grid_subset(grid_ggd, 22)
        @time subset_sol = get_grid_subset(grid_ggd, 23)
        subset_odr = get_grid_subset(grid_ggd, 24)
        subset_idr = get_grid_subset(grid_ggd, 25)
        subset_otarget = get_grid_subset(grid_ggd, 13)
        subset_itarget = get_grid_subset(grid_ggd, 14)

        print("get_space(subset_core) time: ")
        @time space_copy = get_space(subset_core)
        @test space == space_copy
        print("get_space(subset_sol) time: ")
        @time space_copy = get_space(subset_sol)
        @test space == space_copy

        print("get_grid_ggd(subset_core) time: ")
        @time grid_ggd_copy = get_grid_ggd(subset_core)
        @test grid_ggd == grid_ggd_copy
        print("get_grid_ggd(subset_sol) time: ")
        @time grid_ggd_copy = get_grid_ggd(subset_sol)
        @test grid_ggd == grid_ggd_copy

        subset_corebnd = get_grid_subset(grid_ggd, 15)
        subset_separatrix = get_grid_subset(grid_ggd, 16)
        subset_pfrcut = get_grid_subset(grid_ggd, 8)
        subset_otsep = get_grid_subset(grid_ggd, 103)
        subset_itsep = get_grid_subset(grid_ggd, 104)

        print("get_subset_boundary(space, subset_core) time: ")
        @time core_bdry = get_subset_boundary(space, subset_core)
        print("get_subset_boundary(space, subset_sol) time: ")
        @time sol_bdry = get_subset_boundary(space, subset_sol)
        idr_bdry = get_subset_boundary(space, subset_idr)
        odr_bdry = get_subset_boundary(space, subset_odr)

        print("subset_do(intersect, idr_bdry, odr_bdry) time: ")
        @time subset_pfrcut_copy = subset_do(intersect, idr_bdry, odr_bdry)
        @test subset_pfrcut.element == subset_pfrcut_copy.element

        print("subset_do(setdiff, core_bdry, sol_bdry) time: ")
        @time subset_corebnd_copy = subset_do(setdiff, core_bdry, sol_bdry)
        @test subset_corebnd.element == subset_corebnd_copy.element
        @test subset_separatrix.element ==
              subset_do(intersect, sol_bdry,
            subset_do(union, core_bdry, odr_bdry, idr_bdry)).element
        @test subset_otsep.element ==
              subset_do(
            intersect,
            subset_separatrix,
            subset_otarget;
            use_nodes=true,
        ).element
        @test subset_itsep.element ==
              subset_do(
            intersect,
            subset_separatrix,
            subset_itarget;
            use_nodes=true,
        ).element

        print("deepcopy_subset(subset_sol) time: ")
        @time subset_sol_copy = deepcopy_subset(subset_sol)
        @test subset_sol == subset_sol_copy
    end
end

if isempty(ARGS) || "in" in ARGS
    @testset "Test ∈" begin
        grid_ggd = ids.edge_profiles.grid_ggd[1]
        space = grid_ggd.space[1]
        subset_corebnd = get_grid_subset(grid_ggd, 15)
        subset_sol = get_grid_subset(grid_ggd, 23)
        subset_odr = get_grid_subset(grid_ggd, 24)

        @test (6.0, 0.0) ∈ (subset_corebnd, space)
        @test (5.0, -2.5) ∉ (subset_corebnd, space)
        @test (6.0, 4.0) ∈ (subset_sol, space)
        @test (6.0, 3.0) ∉ (subset_sol, space)
        @test (5.1, -3.7) ∈ (subset_odr, space)
        @test (4.5, -3.7) ∉ (subset_odr, space)
    end
end

if isempty(ARGS) || "types" in ARGS
    @testset "Test types" begin
        grid_ggd = ids.edge_profiles.grid_ggd[1]
        resize!(ids.radiation.grid_ggd, 1)
        ids.radiation.grid_ggd[1].path = "edge_profiles/grid_ggd(1)"
        @test grid_ggd.grid_subset == ids.radiation.grid_ggd[1].grid_subset
        @test grid_ggd.identifier == ids.radiation.grid_ggd[1].identifier
        @test grid_ggd.space == ids.radiation.grid_ggd[1].space
        @test grid_ggd.time == ids.radiation.grid_ggd[1].time
    end
end
