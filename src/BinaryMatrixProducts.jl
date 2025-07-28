# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

module BinaryMatrixProducts

export RSMInt, RowSwitchMatrix
include("RowSwitchMatrix.jl")

export ReversibleCircuit, add_gate!, invert_circuit
include("ReversibleCircuit.jl")

export BMP, projbmp, bonddims, volume
export evalfunc
export save_bmp, load_bmp
include("BMP.jl")
export clean1_rl, clean1_lr, clean1
include("BMP_clean.jl")
export apply, minapply, multiapply, layerapply
include("BMP_apply.jl")
include("CustomSparseMat.jl")
include("BMP_minapply.jl")
include("BMP_multiapply.jl")
export insert_var, restrict, erase_var
export compose
export joinfuncs
export swap!, reorder!
include("BMP_ops.jl")

export generate_bmp
include("BMP_util.jl")

export BDD, volume, reduce_bdd!, swap!, save_bdd, load_bdd
include("BDD.jl")
include("Conversions.jl")

export Chip, apply_gate!, minapply_gate!, apply_circuit!, join_chip
include("Chip.jl")

export brute_force!, exact_minimize!, sift!
include("BMP_ordering.jl")
include("BMP_ordering_debug.jl")

export WordBMP
include("WordBMP.jl")

include("TreeStructuredCircuit.jl")

end
