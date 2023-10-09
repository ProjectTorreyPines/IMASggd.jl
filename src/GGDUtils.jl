module GGDUtils

using OMAS: OMAS

export project_prop_on_subset!
export get_subset_centers
export get_prop_with_grid_subset_index

include("types.jl")

include("subset_tools.jl")

include("interpolations.jl")

include("recipes.jl")

end # module GGDUtils
