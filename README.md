# IMASggd

![Format Check](https://github.com/ProjectTorreyPines/IMASggd.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/IMASggd.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/IMASggd.jl/actions/workflows/test.yml/badge.svg)
[![codecov](https://codecov.io/gh/ProjectTorreyPines/IMASggd.jl/graph/badge.svg?token=ZJBRLAXIS1)](https://codecov.io/gh/ProjectTorreyPines/IMASggd.jl)

Package holding utilities for Generalized Grid Description (GGD) objects in IMAS datastructure. It provides, interpolation routines, grid subset tools, and plotting recipes. For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/IMASggd.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/IMASggd.jl/dev)

## Installation

For installation:

```
using Pkg
Pkg.add("IMASggd")
```

## For developers

If you are contributing to this project, you would need to install [dvc](https://dvc.org/) to fetch sample files for testing. Once installed, please configure your ssh so that you can ssh into omega tunneling through cybele without requiring to enter password. This is optional but will make it much easier for you.

Once you have completed above steps, inside the git repo, simply do:
```bash
dvc pull
```

This would download the sample files in the `samples` directory. Then to run tests, you would first need to instantiate the project:
```bash
julia --project
```
Then press `]`:
```julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.11.5 (2025-04-14)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> ]
```
Then type:
```julia
(IMASggd) pkg> instantiate
```
Once the package has been instantiated, you can run the tests using:
```julia
(IMASggd) pkg> test
```
This would run all the tests though. To run specific tests, you can do following from the command line to see help options (this works after you ahve instantiated the project like mentioned above):
```bash
% julia --project test/test.jl help
Usage (from inside IMASggd.jl): 
julia --project test/test.jl [interp] [projection] [in] [subset_tools] [types] [h] [help]

Run tests. Default is all tests.

Optional arguments:
    interp       Test interp                    
    projection   Test project_prop_on_subset!() 
    in           Test ∈                         
    subset_tools Test subset tools              
    types        Test types                     
    h            Show this help message and exit
    help         Show this help message and exit

```
