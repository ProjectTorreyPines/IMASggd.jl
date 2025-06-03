
# IMASggd.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

For installation:

```
using Pkg
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
get_subset_space_objects
get_grid_subset
get_subset_boundary_inds
get_subset_boundary
get_grid_ggd
get_space
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

This function has been used to create following types that are used in this module and
can be imported for further use.

```@docs
IMASggd.all__grid_ggd
IMASggd.all__space
IMASggd.all__grid_subset
IMASggd.all__ggd
IMASggd.all__grid_subset_prop
```

## Plot recipes

Several plot recipes have been defined for easy visualization.
```@docs
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.all__space)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.all__space, ::IMASggd.all__grid_subset)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.all__grid_ggd, ::IMASggd.all__grid_subset_prop)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::AbstractVector{<:IMASggd.all__grid_ggd}, ::IMASggd.all__grid_subset_prop)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer__channel)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer__channel___line_of_sight)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer__channel___n_e_line)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::IMASggd.IMASdd.interferometer__channel___n_e_line_average)
```