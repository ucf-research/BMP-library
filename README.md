# BMP: Binary Matrix Products
BMPs are compact tensor-network representations of Boolean functions. They
are closely related to ordered binary decision diagrams (OBDD). This repository
provides an implementation of basic BMP functionality and is written in Julia.

A detailed explanation of BMPs, their properties, and operations involving BMPs can be found at [https://arxiv.org/abs/2505.01930]

In order to use the code here, you should first create your own Julia environment to
work with:
```
$ mkdir your-project
$ cd your-project
$ julia
julia> ]
pkg> activate .
```
In the same interface, add the package on this repository to your environment:
```
pkg> add https://github.com/ucf-research/BMP-library.git
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

## Documentation
In order to locally compile the project documentation, clone the repository and
run `make.jl` under the `docs` folder:
```
$ git clone https://github.com/ucf-research/BMP-library.git
$ cd BMP-library/docs
$ julia --project make.jl
```
The documentation can then be accessed through `index.html` under
`BMP-library/docs/build`.

## Acknowledgements
This work is based on research supported by the National Science Foundation under grants OIA-2428487 and OIA-2428488.
