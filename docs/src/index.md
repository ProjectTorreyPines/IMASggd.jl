
# IMASggd.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

IMASggd is registered with public repository [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). For installation:

```
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("IMASggd)
```

## Interpolations

Several interpolation functions are available to create interpolaiton functions for data present in a GGD represented over a particular grid subset:

```@docs
interp
get_TPS_mats
get_kdtree
```

## Subset Tools

```@docs
add_subset_element!
get_subset_space
get_grid_subset
get_subset_boundary_inds
get_subset_boundary
subset_do
get_subset_centers
project_prop_on_subset!
deepcopy_subset
Base.:∈
get_prop_with_grid_subset_index
```

## Types

```@docs
get_types_with
```

## Plot recipes

Several plot recipes have been defined for easy visualization.
```@docs
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.edge_profiles__grid_ggd___space)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.edge_profiles__grid_ggd___space, ::IMASggd.IMASdd.edge_profiles__grid_ggd___grid_subset)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.edge_profiles__grid_ggd, ::IMASggd.IMASdd.IDSvectorElement)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::AbstractVector{<:IMASggd.IMASdd.edge_profiles__grid_ggd}, ::IMASggd.IMASdd.IDSvectorElement)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer__channel)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer__channel___line_of_sight)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer__channel___n_e_line)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer__channel___n_e_line_average)
```