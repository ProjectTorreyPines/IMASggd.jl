# using Pkg
# Pkg.add("Documenter")
# Pkg.add(; url="https://github.com/ProjectTorreyPines/IMASdd.jl.git");
# Pkg.develop(PackageSpec(; path="../"));
# Pkg.instantiate()
using Documenter
using IMASggd
using RecipesBase: RecipesBase

makedocs(;
    modules=[IMASggd],
    format=Documenter.HTML(),
    sitename="IMASggd",
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ProjectTorreyPines/IMASggd.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "v#.#"],
)
