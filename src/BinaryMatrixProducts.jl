module BinaryMatrixProducts

export RSMInt, RowSwitchMatrix
include("RowSwitchMatrix.jl")

export Circuit, add_gate!
include("Circuit.jl")

export BMP, projbmp, bonddims, volume
export evalfunc
export clean1_rl, clean1_lr, clean1
export apply, minapply, multiapply, layerapply
export joinfuncs, swap!
include("BMP.jl")

export generate_bmp
include("BMP_util.jl")

export BDD, volume, reduce_bdd!, swap!, save_bdd, load_bdd
include("BDD.jl")
include("Conversions.jl")

export Chip, apply_gate!, minapply_gate!, join_chip
include("Chip.jl")

export exact_minimize!, sift!
include("BMP_ordering.jl")

export WordBMP
include("WordBMP.jl")

include("TreeStructuredCircuit.jl")

end
