# GGDUtils

![Format Check](https://github.com/ProjectTorreyPines/GGDUtils.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/GGDUtils.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/GGDUtils.jl/actions/workflows/test.yml/badge.svg)

Package holding utilities for Generalized Grid Description (GGD) objects in IMAS datastructure. Primary goals are interpolation and core profile extrapolation. For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/GGDUtils.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/GGDUtils.jl/dev)

## Installation

GGDUtils is registered with public repository [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). For installation:

```
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("GGDUtils")
```
