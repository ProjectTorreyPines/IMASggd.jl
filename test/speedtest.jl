import IMASggd: interp, get_kdtree, project_prop_on_subset!, get_grid_subset, IMASdd

println("-----------------------------------------------------------------------------")
print("json2imas() time with compilation: ")
@time ids = IMASdd.json2imas(
    "$(@__DIR__)/../samples/time_dep_edge_profiles_last_step_only.json",
)
print("json2imas() time (true runtime): ")
@time ids = IMASdd.json2imas(
    "$(@__DIR__)/../samples/time_dep_edge_profiles_last_step_only.json",
)

println("-----------------------------------------------------------------------------")
n_e = ids.edge_profiles.ggd[1].electrons.density[1];
grid_ggd = ids.edge_profiles.grid_ggd[1];
space = grid_ggd.space[1]
print("interp(prop, grid_ggd) time with compilation: ")
@time get_n_e = interp(n_e, grid_ggd)
print("interp(prop, grid_ggd) (true runtime): ")
@time get_n_e = interp(n_e, grid_ggd)

println("-----------------------------------------------------------------------------")
subset = get_grid_subset(grid_ggd, -5)

print("interp(prop_arr, space, subset) time with compilation: ")
@time get_n_e = interp(ids.edge_profiles.ggd[1].electrons.density, space, subset)
print("interp(prop_arr, space, subset) time (true runtime): ")
@time get_n_e = interp(ids.edge_profiles.ggd[1].electrons.density, space, subset)

println("-----------------------------------------------------------------------------")
print("interp(prop_arr, grid_ggd, grid_subset_index) time with compilation: ")
@time get_n_e = interp(ids.edge_profiles.ggd[1].electrons.density, grid_ggd, -5)
print("interp(prop_arr, grid_ggd, grid_subset_index) time (true runtime): ")
@time get_n_e = interp(ids.edge_profiles.ggd[1].electrons.density, grid_ggd, -5)

println("-----------------------------------------------------------------------------")
kdtree = get_kdtree(space)
print("interp(prop_arr, kdtree) time with compilation: ")
@time get_T_e = interp(ids.edge_profiles.ggd[1].electrons.temperature[1].values, kdtree)
print("interp(prop_arr, kdtree) time (true runtime): ")
@time get_T_e = interp(ids.edge_profiles.ggd[1].electrons.temperature[1].values, kdtree)

println("-----------------------------------------------------------------------------")
prop = ids.edge_profiles.ggd[1].electrons.density
# All cells
from_subset = get_grid_subset(ids.edge_profiles.grid_ggd[1], -5)
# separatix
to_subset = get_grid_subset(ids.edge_profiles.grid_ggd[1], 16)
print(
    "project_prop_on_subset!(prop, from_subset, to_subset, space) time with compilation: ",
)
@time separatix_centers, values_at_separatix =
    project_prop_on_subset!(prop, from_subset, to_subset, space)
print(
    "project_prop_on_subset!(prop, from_subset, to_subset, space) time (true runtime): ",
)
@time separatix_centers, values_at_separatix =
    project_prop_on_subset!(prop, from_subset, to_subset, space)

println("-----------------------------------------------------------------------------")
subset_core =
    get_grid_subset(ids.edge_profiles.grid_ggd[1], 22)
print("project_prop_on_subset!(prop, from_subset, subset_core) time with compilation: ")
@time core_element_inds, values_at_core =
    project_prop_on_subset!(prop, from_subset, subset_core)
print("project_prop_on_subset!(prop, from_subset, subset_core) time (true runtime): ")
@time core_element_inds, values_at_core =
    project_prop_on_subset!(prop, from_subset, subset_core)

println("-----------------------------------------------------------------------------")
subset_corebnd = get_grid_subset(grid_ggd, 15)
subset_sol = get_grid_subset(grid_ggd, 23)
subset_odr = get_grid_subset(grid_ggd, 24)

print("test ∈ (edges) time with compilation: ")
@time (6.0, 0.0) ∈ (subset_corebnd, space)
print("test ∈ (edges) time (true runtime): ")
@time (6.0, 0.0) ∈ (subset_corebnd, space)

print("test ∈ (cells) time with compilation: ")
@time (6.0, 4.0) ∈ (subset_sol, space)
print("test ∈ (cells) time (true runtime): ")
@time (6.0, 4.0) ∈ (subset_sol, space)
