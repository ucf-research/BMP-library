# BMP: Binary Matrix Products
BMPs are compact representations of boolean functions as tensor networks. They
are closely related to ordered binary decision diagrams (OBDD). This repository
provides an implementation of basic BMP functionality written in Julia.

In order to use the code here, you must first clone the repository:
```
$ git clone https://github.com/ucf-research/BMP-library.git
```
Once you have a local copy of the package, create your own Julia environment to
work with:
```
$ mkdir your-project
$ cd your-project
$ julia
julia> ]
pkg> activate .
```
In the same interface, add the local copy of the package to your environment:
```
pkg> add /path/to/BMP-library
```
Now you can import the BMP functions by adding the following line to your Julia
files:
```
using BinaryMatrixProducts
```
In order to run these files you need to use the environment you created:
```
$ julia --project=/path/to/your-project foo.jl
```

This work is based on research supported by the National Science Foundation under grants OIA-2428487 and OIA-2428488.
