# BMP Library Documentation homepage

This is the homepage of the documentation for the BMP library. If you want to
use the library in your Julia code, you need to add it to your project as
follows:
```
$ mkdir your-project
$ cd your-project
$ julia --project
julia> using Pkg
julia> Pkg.add(url="https://github.com/ucf-research/BMP-library.git")
```
Your Julia scripts using the library should include the line
```julia
using BinaryMatrixProducts
```
You can run these scripts using the command below.
```
$ julia --project=/path/to/your-project foo.jl
```
