# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

"""
    Chip

A data type that represents the BMP of a reversible function ``\\{0,1\\}^n \\to \\{0,1\\}^n``.
Each of the `n` bitlines is stored as a `BareBMP` that represents the output on that bitline
as a function of all the primary inputs. These BMPs are constrained to have the same variable
ordering and to be in canonical form at all times for easier processing.
"""
struct Chip
    bitlines::Vector{BareBMP}
    order::Vector{BMPVarInt}
    position::Vector{BMPVarInt}
    function Chip(n::Integer, order)
        position = Vector{BMPVarInt}(undef, length(order))
        position[order] .= 1:length(order)
        bitlines = [projbmp_bare(xi, n) for xi in position]
        return new(bitlines, order, position)
    end
end

"""
    Chip(n::Integer)

Creates a `Chip` for the identity permutation on `n` bitlines.
"""
function Chip(n::Integer)
    return Chip(n, collect(1:n))
end

"""
    Chip(circuit::ReversibleCircuit)

Converts `circuit` to a `Chip`. This may be a costly operation
depending on the circuit.
"""
function Chip(circuit::ReversibleCircuit)
    chip = Chip(circuit.n, collect(1:circuit.n))
    for cg in circuit.gates
        minapply_gate!(chip, cg)
    end
    return chip
end

"""
    volume(chip::Chip)

Returns the sum of the volumes of the BMPs associated with the bitlines.
"""
function volume(chip::Chip)
    cnt = 0
    for mats in chip.bitlines
        cnt += sum(length(mats[i,1].rows) for i in axes(mats,1))
    end
    return cnt
end

"""
    evalfunc(chip::Chip, input::AbstractArray)

Computes the output of the permutation represented by `chip` on inputs given
in `input`.
"""
function evalfunc(chip::Chip, input::AbstractArray)
    n = size(input, 1)
    n_samps = div(length(input), size(input, 1))
    result = similar(input, (n, n_samps))
    for (i,bl) in enumerate(chip.bitlines)
        temp = evalfunc(bl, input, [0,1], chip.order)
        temp = reshape(temp, (1, n_samps))
        result[i,:] .= temp[1,:]
    end
    return reshape(result, size(input))
end

"""
    apply_gate!(chip::Chip, tab::Array{<:Integer}, bits)

Realizes the action of a reversible gate given by `tab` on bitlines `bits` of
`chip`. This method uses the direct-product APPLY operation internally.
"""
function apply_gate!(chip::Chip, tab::Array{<:Integer}, bits)
    tensor_bmp = apply_mats(map(i -> chip.bitlines[i], bits))
    n_bits = length(bits)
    for (i,b) in zip(n_bits-1:-1:0, bits)
        bit_tab = tab .>> i .& 1
        chip.bitlines[b] = clean1(tensor_bmp, bit_tab)
    end
end

"""
    apply_gate!(chip::Chip, g::ReversibleGate)

Realizes the action of gate `g` on `chip`, using the direct-product APPLY operation
internally.
"""
function apply_gate!(chip::Chip, g::ReversibleGate)
    apply_gate!(chip, g.perm, g.bits)
end

"""
    apply_circuit!(chip::Chip, circuit::ReversibleCircuit)

Realizes the action of all the gates in `circuit` on `chip`. This method uses
the direct-product APPLY operation internally.
"""
function apply_circuit!(chip::Chip, circuit::ReversibleCircuit)
    if circuit.n != length(chip.bitlines)
        throw(DimensionMismatch("Number of bitlines on the chip and the circuit don't match."))
    end
    for cg in circuit.gates
        apply_gate!(chip, cg)
    end
end

"""
    minapply_gate!(chip::Chip, tab::Array{<:Integer}, bits)

Realizes the action of a reversible gate given by `tab` on bitlines `bits` of
`chip`. This method uses the direct-sum APPLY operation internally.
"""
function minapply_gate!(chip::Chip, tab::Array{<:Integer}, bits)
    sum_bmp, U = minapply_mats(map(i -> chip.bitlines[i], bits))
    n_bits = length(bits)
    Rs = map(i -> [0,1], bits)
    for (i,b) in zip(n_bits-1:-1:0, bits)
        bit_tab = tab .>> i .& 1
        R = minapply_term(bit_tab, U, Rs)
        chip.bitlines[b] = clean1_rl(sum_bmp, R)
    end
end

"""
    minapply_gate!(chip::Chip, g::ReversibleGate)

Realizes the action of gate `g` on `chip`, using the direct-sum APPLY operation
internally.
"""
function minapply_gate!(chip::Chip, g::ReversibleGate)
    minapply_gate!(chip, g.perm, g.bits)
end

"""
    minapply_circuit!(chip::Chip, circuit::ReversibleCircuit)

Realizes the action of all the gates in `circuit` on `chip`. This method uses
the direct-sum APPLY operation internally.
"""
function minapply_circuit!(chip::Chip, circuit::ReversibleCircuit)
    if circuit.n != length(chip.bitlines)
        throw(DimensionMismatch("Number of bitlines on the chip and the circuit don't match."))
    end
    for cg in circuit.gates
        minapply_gate!(chip, cg)
    end
end

"""
    join_chip(chip::Chip)

Combines all the bitlines of `chip` into a standard BMP with multiple outputs.
The BMP inherits its input variable ordering from `chip.` On the other hand,
the outputs are ordered according to their labels.
"""
function join_chip(chip::Chip)
    n = length(chip.bitlines)
    mats = Matrix{RowSwitchMatrix}(undef, (n, 2))
    for i=1:n, j=1:2
        mats[i,j] = dsum([bmp[i,j] for bmp in chip.bitlines])
    end
    R = fill(0, 2*n)
    R[2:2:end] .= 1
    return clean1(BMP(mats, R, copy(chip.order)))
end
