
# GGDUtils.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

### Using make:
After cloning this repo, check the make menu:
```
GGDUtils.jl % make help
Help Menu

make env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories
make env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones
make clean: Deletes Project.toml and Manifest.toml for a fresh start
```

#### make r
This option creates local copies of required private repositories at the same level as current repository and uses them in develop mode to create a Manifest.toml

#### make u
This option uses url of required private repositories to create a static Manifest.toml attached to current master branches of these repositories.

#### make clean
Deletes Manifest.toml so that environment can be recreated, to update or change the last used method.

### Using Julia REPL and installing using Github url

Or, in julia REPL:
```julia
julia> using Pkg;
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/IMASDD.jl.git");
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/GGDUtils.jl.git");
julia> Pkg.instantiate()
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
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.edge_profiles__grid_ggd___space)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.edge_profiles__grid_ggd___space, ::GGDUtils.IMASDD.edge_profiles__grid_ggd___grid_subset)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.edge_profiles__grid_ggd, ::GGDUtils.IMASDD.IDSvectorElement)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::AbstractVector{<:GGDUtils.IMASDD.edge_profiles__grid_ggd}, ::GGDUtils.IMASDD.IDSvectorElement)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer__channel)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer__channel___line_of_sight)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer__channel___n_e_line)
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer__channel___n_e_line_average)
```