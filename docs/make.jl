# using Pkg
# Pkg.add("Documenter")
# Pkg.add(; url="https://github.com/ProjectTorreyPines/IMASdd.jl.git");
# Pkg.develop(PackageSpec(; path="../"));
# Pkg.instantiate()
using Documenter
using GGDUtils
using RecipesBase: RecipesBase

makedocs(;
    modules=[GGDUtils],
    format=Documenter.HTML(),
    sitename="GGDUtils",
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ProjectTorreyPines/GGDUtils.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "v#.#"],
)
