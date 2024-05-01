var documenterSearchIndex = {"docs":
[{"location":"#GGDUtils.jl","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"Pages = [\"index.md\"]\nDepth = 5","category":"page"},{"location":"#Installation","page":"GGDUtils.jl","title":"Installation","text":"","category":"section"},{"location":"#Using-make:","page":"GGDUtils.jl","title":"Using make:","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"After cloning this repo, check the make menu:","category":"page"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"GGDUtils.jl % make help\nHelp Menu\n\nmake env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories\nmake env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones\nmake clean: Deletes Project.toml and Manifest.toml for a fresh start","category":"page"},{"location":"#make-r","page":"GGDUtils.jl","title":"make r","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"This option creates local copies of required private repositories at the same level as current repository and uses them in develop mode to create a Manifest.toml","category":"page"},{"location":"#make-u","page":"GGDUtils.jl","title":"make u","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"This option uses url of required private repositories to create a static Manifest.toml attached to current master branches of these repositories.","category":"page"},{"location":"#make-clean","page":"GGDUtils.jl","title":"make clean","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"Deletes Manifest.toml so that environment can be recreated, to update or change the last used method.","category":"page"},{"location":"#Using-Julia-REPL-and-installing-using-Github-url","page":"GGDUtils.jl","title":"Using Julia REPL and installing using Github url","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"Or, in julia REPL:","category":"page"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"julia> using Pkg;\njulia> Pkg.add(; url=\"https://github.com/ProjectTorreyPines/IMASDD.jl.git\");\njulia> Pkg.add(; url=\"https://github.com/ProjectTorreyPines/GGDUtils.jl.git\");\njulia> Pkg.instantiate()","category":"page"},{"location":"#Interpolations","page":"GGDUtils.jl","title":"Interpolations","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"Several interpolation functions are available to create interpolaiton functions for data present in a GGD represented over a particular grid subset:","category":"page"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"interp\nget_TPS_mats\nget_kdtree","category":"page"},{"location":"#GGDUtils.interp","page":"GGDUtils.jl","title":"GGDUtils.interp","text":"interp(\nprop_values::Vector{T},\nkdtree::KDTree;\nuse_nearest_n::Int=4,\nweighing::Function=(d) -> 1 / d,\n\n) where {T <: Real}\n\nLowest level interpolation function. It takes a vector of property values and a KDTree defined over a 2D space with the same number of nodes as the property values. It returns a function that can be used to interpolate the property values at any point in the space.\n\n\n\n\n\ninterp(\ny::Vector{T},\nTPS_mats::Tuple{Matrix{U}, Matrix{U}, Matrix{U}, Vector{Tuple{U, U}}},\n\n) where {T <: Real, U <: Real}\n\nLowest level function for Thin Plate Spline method\n\n\n\n\n\ninterp(y::Vector{T}, x::Vector{Tuple{U, U}}) where {T <: Real, U <: Real}\n\nThin plate smoothing interpolation function for a 2d space scalar function. The algorithm has been adopted from:\n\nhttp://www.geometrictools.com/Documentation/ThinPlateSplines.pdf\n\nThis is an implementation of Euler-Lagrange equation for minimizing bending energy of a surface.\n\n\n\n\n\ninterp(\nprop_values::Vector{T},\nspace::IMASDD.edge_profiles__grid_ggd___space\n\n) where {T <: Real}\n\nIf the whole space is provided instead of a kdtree, calculate the kdtree for whole space. Again, here it is assumed that the property values are porvided for each node of the space.\n\n\n\n\n\ninterp(\nprop_values::Vector{Real},\nspace::IMASDD.edge_profiles__grid_ggd___space,\nsubset::IMASDD.edge_profiles__grid_ggd___grid_subset\n\n)\n\nIf a subset of the space is provided, calculate the kdtree for the subset. In this case it is assumed that the property values are provided for each element of the subset.\n\n\n\n\n\ninterp(\nprop::edge_profiles__prop_on_subset,\ngrid_ggd::IMASDD.edge_profiles__grid_ggd,\nvalue_field::Symbol=:values\n\n)\n\nExample: gridggd = dd.edgeprofiles.gridggd[1] getelectrondensity = interp(dd.edgeprofiles.ggd[1].electrons.density[1], gridggd) getefieldpar = interp(dd.edgeprofiles.ggd[1].efield[1], grid_ggd, :parallel)\n\n\n\n\n\ninterp(\nprop_arr::AbstractVector{T},\nspace::IMASDD.edge_profiles__grid_ggd___space,\nsubset::IMASDD.edge_profiles__grid_ggd___grid_subset,\nvalue_field::Symbol=:values\n\n) where {T <: edgeprofiles__propon_subset}\n\nExample: sol = getgridsubsetwithindex(dd.edgeprofiles.gridggd[1], 23) getelectrondensity = interp(dd.edge_profiles.ggd[1].electrons.density, space, sol)\n\n\n\n\n\ninterp(\nprop_arr::AbstractVector{T},\ngrid_ggd::IMASDD.edge_profiles__grid_ggd,\ngrid_subset_index::Int,\nvalue_field::Symbol=:values\n\n) where {T <: edgeprofiles__propon_subset}\n\nExample: getnesep = interp(dd.edgeprofiles.ggd[1].electrons.density, grid_ggd, 16)\n\n\n\n\n\ninterp(eqt::IMASDD.equilibrium__time_slice)\n\nFor a given equilibrium time slice, return a function that can be used to interpolate from (r, z) space to rho (normalized toroidal flux coordinate)space.\n\nExample: rz2rho = interp(dd.equilibrium.time_slice[1]) rho = rz2rho.([(r, z) for r in 3:0.01:9, for z in -5:0.01:5])\n\n\n\n\n\ninterp(\nprop::Vector{T},\nprof::IMASDD.core_profiles__profiles_1d,\n\n) where {T <: Real}\n\nReturns an inteprolation function for the core profile property values defined on normalized toroidal flux coordinate rho.\n\nExample: coreprofilene = dd.coreprofiles.profiles1d[1].electrons.density getne = interp(coreprofilene, dd.coreprofiles.profiles1d[1]) getne(1) # Returns electron density at rho = 1 (separatix)\n\n\n\n\n\ninterp(\nprop::Vector{T},\nprof::IMASDD.core_profiles__profiles_1d,\nrz2rho::Function,\n\n)\n\nReturns an inteprolation function in (R, Z) domain for the core profile property values defined on normalized toroidal flux coordinate rho and with a provided function to convert (R,Z) to rho.\n\nExample:\n\nrz2rho = interp(dd.equilibrium.timeslice[1]) coreprofilene = dd.coreprofiles.profiles1d[1].electrons.density getne = interp(coreprofilene, dd.coreprofiles.profiles1d[1], rz2rho) getn_e(5.0, 3.5) # Returns electron density at (R, Z) = (5.0, 3.5)\n\n\n\n\n\ninterp(\nprop::Vector{T},\nprof::IMASDD.core_profiles__profiles_1d,\neqt::IMASDD.equilibrium__time_slice,\n\n) where {T <: Real}\n\nReturns an inteprolation function in (R, Z) domain for the core profile property values defined on normalized toroidal flux coordinate rho and with a provided equilibrium time slice to get (R, Z) to rho conversion.\n\nExample:\n\neqt = dd.equilibrium.timeslice[1] coreprofilene = dd.coreprofiles.profiles1d[1].electrons.density getne = interp(coreprofilene, dd.coreprofiles.profiles1d[1], eqt) getn_e(5.0, 3.5) # Returns electron density at (R, Z) = (5.0, 3.5)\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.get_TPS_mats","page":"GGDUtils.jl","title":"GGDUtils.get_TPS_mats","text":"get_TPS_mats(x::Vector{Tuple{U, U}}) where {U <: Real}\n\nWe can separate matrix operations of thin plate spline method to utilize the calculation for itnerpolating over same space. Similar to the idea of reusing KDTree like above.\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.get_kdtree","page":"GGDUtils.jl","title":"GGDUtils.get_kdtree","text":"get_kdtree(space::IMASDD.edge_profiles__grid_ggd___space)\n\nGet a KDTree for all the cells in the space for search for nearest neighbours.\n\n\n\n\n\nget_kdtree(\nspace::IMASDD.edge_profiles__grid_ggd___space,\nsubset::IMASDD.edge_profiles__grid_ggd___grid_subset,\n\n)\n\nGet a KDTree for a subset of the space for search for nearest neighbours.\n\n\n\n\n\n","category":"function"},{"location":"#Subset-Tools","page":"GGDUtils.jl","title":"Subset Tools","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"add_subset_element!\nget_subset_space\nget_grid_subset\nget_subset_boundary_inds\nget_subset_boundary\nsubset_do\nget_subset_centers\nproject_prop_on_subset!\ndeepcopy_subset\nBase.:∈\nget_prop_with_grid_subset_index","category":"page"},{"location":"#GGDUtils.add_subset_element!","page":"GGDUtils.jl","title":"GGDUtils.add_subset_element!","text":"add_subset_element!(\nsubset::IMASDD.edge_profiles__grid_ggd___grid_subset,\nsn::Int,\ndim::Int,\nindex::Int,\nin_subset=(x...) -> true;\nkwargs...,\n\n)\n\nAdds a new element to girdsubset with properties space number (sn), dimension (dim), and index (index). The element is added only if the function insubset returns true.\n\n\n\n\n\nadd_subset_element!(\nsubset,\nsn,\ndim,\nindex::Vector{Int},\nin_subset=(x...) -> true;\nkwargs...,\n\n)\n\nOverloaded to work differently (faster) with list of indices to be added.\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.get_subset_space","page":"GGDUtils.jl","title":"GGDUtils.get_subset_space","text":"get_subset_space(space::IMASDD.edge_profiles__grid_ggd___space,\nelements::AbstractVector{<:IMASDD.edge_profiles__grid_ggd___grid_subset___element})\n\nReturns an array of space object indices corresponding to the correct objectsperdimension (nodes, edges or cells) for the subset elements.\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.get_grid_subset","page":"GGDUtils.jl","title":"GGDUtils.get_grid_subset","text":"get_grid_subset(\ngrid_ggd::IMASDD.edge_profiles__grid_ggd,\ngrid_subset_index::Int,\n\n)\n\nReturns the gridsubset in a gridggd with the matching gridsubsetindex\n\n\n\n\n\nget_grid_subset(\ngrid_ggd::IMASDD.edge_profiles__grid_ggd,\ngrid_subset_name::String,\n\n)\n\nReturns the gridsubset in a gridggd with the matching gridsubsetname\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.get_subset_boundary_inds","page":"GGDUtils.jl","title":"GGDUtils.get_subset_boundary_inds","text":"get_subset_boundary_inds(\nspace::IMASDD.edge_profiles__grid_ggd___space,\nsubset::IMASDD.edge_profiles__grid_ggd___grid_subset,\n\n)\n\nReturns an array of space object indices corresponding to the boundary of the subset. That means, it returns indices of nodes that are at the end of open edge subset or it returns the indices of edges that are the the boundary of a cell subset. Returns an empty array if the subset is 1D (nodes).\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.get_subset_boundary","page":"GGDUtils.jl","title":"GGDUtils.get_subset_boundary","text":"get_subset_boundary(\nspace::IMASDD.edge_profiles__grid_ggd___space,\nsubset::IMASDD.edge_profiles__grid_ggd___grid_subset,\n\n)\n\nReturns an array of elements of grid_subset generated from the boundary of the subset provided. The dimension of these elments is reduced by 1.\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.subset_do","page":"GGDUtils.jl","title":"GGDUtils.subset_do","text":"subset_do(set_operator,\nitrs::Vararg{AbstractVector{<:IMASDD.edge_profiles__grid_ggd___grid_subset___element}};\nspace::IMASDD.edge_profiles__grid_ggd___space=IMASDD.edge_profiles__grid_ggd___space(),\nuse_nodes=false)\n\nFunction to perform any set operation (intersect, union, setdiff etc.) on subset.element to generate a list of elements to go to subset object. If usenodes is true, the set operation will be applied on the set of nodes from subset.element, space argument is required for this. Note: that the arguments are subset.element (not the subset itself). Similarly, the return object is a list of IMASDD.edgeprofiles__gridggd___gridsubset___element.\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.get_subset_centers","page":"GGDUtils.jl","title":"GGDUtils.get_subset_centers","text":"get_subset_centers(space::IMASDD.edge_profiles__grid_ggd___space,\n                        subset::IMASDD.edge_profiles__grid_ggd___grid_subset)\n\nReturns an array of tuples corresponding to (r,z) coordinates of the center of cells or the center of edges in the subset space.\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.project_prop_on_subset!","page":"GGDUtils.jl","title":"GGDUtils.project_prop_on_subset!","text":"project_prop_on_subset!(prop_arr::AbstractVector{T},\nfrom_subset::IMASDD.edge_profiles__grid_ggd___grid_subset,\nto_subset::IMASDD.edge_profiles__grid_ggd___grid_subset,\nvalue_field::Symbol=:values,\n\n) where {T <: edgeprofiles__propon_subset}\n\nIf the dimensions of fromsubset and tosubset are the same, this function can be used to add another instance on a property vector representing the value in tosubset without any interpolation or use of space object. The function returns a tuple of indices of elements of tosubset and the values of the property in to_subset.\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.deepcopy_subset","page":"GGDUtils.jl","title":"GGDUtils.deepcopy_subset","text":"deepcopy_subset(subset::IMASDD.edge_profiles__grid_ggd___grid_subset)\n\nFaster deepcopy function for gridsubset object. This function is used to create a deep copy of a gridsubset object bypassing several checks performed by IMASDD.\n\n\n\n\n\n","category":"function"},{"location":"#Base.:∈","page":"GGDUtils.jl","title":"Base.:∈","text":"Base.:∈(\npoint::Tuple{Real, Real},\nsubset_of_space::Tuple{\n    IMASDD.edge_profiles__grid_ggd___grid_subset,\n    IMASDD.edge_profiles__grid_ggd___space,\n},\n\n)\n\nOverloading ∈ operator to check if a point is inside a subset of space.\n\nIf the subset is 1-dimensional, all points are searched. If the subset is 2-dimensional, it is checked if the point is within the enclosed area. It is assumed that a 2-dimensional subset used in such a context will form a closed area. If the subset is 3-dimensional, its boundary is calculated on the fly. If used multiple times, it is recommended to calculate the boundary once and store it in a variable.\n\n\n\n\n\n","category":"function"},{"location":"#GGDUtils.get_prop_with_grid_subset_index","page":"GGDUtils.jl","title":"GGDUtils.get_prop_with_grid_subset_index","text":"get_prop_with_grid_subset_index(\nprop::AbstractVector{T},\ngrid_subset_index::Int,\n\n) where {T <: edgeprofiles__propon_subset}\n\nFind the edgeprofiles property instance in an array of properties that corresponds to the gridsubset_index provided.\n\n\n\n\n\n","category":"function"},{"location":"#Types","page":"GGDUtils.jl","title":"Types","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"get_types_with","category":"page"},{"location":"#GGDUtils.get_types_with","page":"GGDUtils.jl","title":"GGDUtils.get_types_with","text":"get_types_with(parent::Type, field::Symbol)\n\nA type creation utility meant for searching types in IMAS database. This function returns a list of types that are fields at any level below the parent data type which have a particular field present in it.\n\nExample:\n\ngettypeswith(IMASDD.edgeprofiles, :gridsubset_index)\n\nreturns all edgeprofiles types that have a subfield named gridsubset_index.\n\n\n\n\n\n","category":"function"},{"location":"#Plot-recipes","page":"GGDUtils.jl","title":"Plot recipes","text":"","category":"section"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"Several plot recipes have been defined for easy visualization.","category":"page"},{"location":"","page":"GGDUtils.jl","title":"GGDUtils.jl","text":"RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.edge_profiles__grid_ggd___space)\nRecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.edge_profiles__grid_ggd___space, ::GGDUtils.IMASDD.edge_profiles__grid_ggd___grid_subset)\nRecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.edge_profiles__grid_ggd, ::GGDUtils.IMASDD.IDSvectorElement)\nRecipesBase.apply_recipe(::Dict{Symbol,Any}, ::AbstractVector{<:GGDUtils.IMASDD.edge_profiles__grid_ggd}, ::GGDUtils.IMASDD.IDSvectorElement)\nRecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer)\nRecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer__channel)\nRecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer__channel___line_of_sight)\nRecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer__channel___n_e_line)\nRecipesBase.apply_recipe(::Dict{Symbol,Any}, ::GGDUtils.IMASDD.interferometer__channel___n_e_line_average)","category":"page"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, IMASDD.edge_profiles__grid_ggd___space}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(space::IMASDD.edge_profiles__grid_ggd___space)\n\nPlot the grid_ggd space object. Defaults to size of [600, 900] and linecolor of :black, linewidth of 0.2, and no legend.\n\n\n\n\n\n","category":"method"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, IMASDD.edge_profiles__grid_ggd___space, IMASDD.edge_profiles__grid_ggd___grid_subset}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(\nspace::IMASDD.edge_profiles__grid_ggd___space,\nsubset::IMASDD.edge_profiles__grid_ggd___grid_subset,\n\n)\n\nPlot the a subset of a space. Defaults to size of [600, 900] and linecolor of :black, linewidth of 0.2, and no legend.\n\n\n\n\n\n","category":"method"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, IMASDD.edge_profiles__grid_ggd, IMASDD.IDSvectorElement}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(\ngrid_ggd::IMASDD.edge_profiles__grid_ggd,\nprop::IMASDD.IDSvectorElement,\n\n)\n\nPlot 2D heatmap of edgeprofilesggd property on a gridggd space object. Defaults to size of [635, 900], xaxis of \"R / m\", yaxis of \"Z / m\", and no legend. If :seriescolor is not provided, :inferno color scheme is used. If :colorbartitle is not provided, the property name is used. This function creates a plot with layout [a{0.95w} b] where a is the heatmap and b is the colorbar.\n\n\n\n\n\n","category":"method"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, AbstractVector{<:IMASDD.edge_profiles__grid_ggd}, IMASDD.IDSvectorElement}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(\ngrid_ggd_arr::AbstractVector{<:IMASDD.edge_profiles__grid_ggd},\nprop::IMASDD.IDSvectorElement,\n\n)\n\nPlot 2D heatmap of edgeprofilesggd property on a gridggd space object. Defaults to size of [635, 900], xaxis of \"R / m\", yaxis of \"Z / m\", and no legend. If :seriescolor is not provided, :inferno color scheme is used. If :colorbartitle is not provided, the property name is used. This function creates a plot with layout [a{0.95w} b] where a is the heatmap and b is the colorbar.\n\n\n\n\n\n","category":"method"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, IMASDD.interferometer}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(\nifo::IMASDD.interferometer,\n\n)\n\nPlot all the channels of interferometer object. Optional keywords:\n\n:plottype: :los(default), :ne, or :neaverage. :los plots the line of sight of the channel in a 2D plot, :ne plots the integrated ne along the line of sight vs time, and :neaverage plots the average n_e vs time.\n:mirror: true(default) or false.\n:mirror_length: 0.5(default).\n:mirror_thickness: 0.1(default).\n\n\n\n\n\n","category":"method"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, IMASDD.interferometer__channel}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(\nifo_ch::IMASDD.interferometer__channel,\n\n)\n\nPlot individual channel of interferometer. Optional keywords:\n\n:plottype: :los(default), :ne, or :neaverage. :los plots the line of sight of the channel in a 2D plot, :ne plots the integrated ne along the line of sight vs time, and :neaverage plots the average n_e vs time.\n:mirror: true(default) or false.\n:mirror_length: 0.5(default).\n:mirror_thickness: 0.1(default).\n\n\n\n\n\n","category":"method"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, IMASDD.interferometer__channel___line_of_sight}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(\nifo_ch_los::IMASDD.interferometer__channel___line_of_sight,\n\n)\n\nPlot line of sight of a channel of interferometer.\n\nDefault plot settings:\n\nsubplot: 1\nsize: [600, 900]\nxaxis: \"R / m\"\nyaxis: \"Z / m\"\n\nOptional keywords:\n\n:mirror: true(default) or false.\n:mirror_length: 0.5(default).\n:mirror_thickness: 0.1(default).\n\n\n\n\n\n","category":"method"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, IMASDD.interferometer__channel___n_e_line}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(\nifo_ch_n_e_line::IMASDD.interferometer__channel___n_e_line,\n\n)\n\nPlot line integrated electron density of a channel of interferometer.\n\nDefault plot settings:\n\nsubplot: 1\nxaxis: \"time / s\"\nyaxis: \"Integerated n_e / m^-2\"\n\nOptional keywords:\n\n:average: true or false(default). If true, plot the average n_e vs time.\n\n\n\n\n\n","category":"method"},{"location":"#RecipesBase.apply_recipe-Tuple{Dict{Symbol, Any}, IMASDD.interferometer__channel___n_e_line_average}","page":"GGDUtils.jl","title":"RecipesBase.apply_recipe","text":"plot(\nifo_ch_n_e_line_average::IMASDD.interferometer__channel___n_e_line_average,\n\n)\n\nPlot average electron density of a channel of interferometer.\n\nDefault plot settings:\n\nsubplot: 1\nxaxis: \"time / s\"\nyaxis: \"Average n_e / m^-3\"\n\n\n\n\n\n","category":"method"}]
}