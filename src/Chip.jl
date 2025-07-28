# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

struct Chip
    bitlines::Vector{BareBMP}
    order::Vector{UInt32}
    position::Vector{UInt32}
    function Chip(n::Integer, order::Vector{<:Integer})
        position = fill(UInt32(0), length(order))
        position[order] .= order
        bitlines = [projbmp_bare(xi, n) for xi in position]
        return new(bitlines, order, position)
    end
end

function Chip(n::Integer)
    return Chip(n, collect(1:n))
end

function Chip(circuit::ReversibleCircuit)
    chip = Chip(circuit.n, collect(1:circuit.n))
    for cg in circuit.gates
        apply_gate!(chip, cg)
    end
    return chip
end

function evalfunc(chip::Chip, input::BitArray)
    n = size(input, 1)
    n_samps = div(length(input), size(input, 1))
    result = BitArray(undef, (n, n_samps))
    for (i,bl) in enumerate(chip.bitlines)
        temp = evalfunc(bl, input, [0,1], chip.order)
        temp = reshape(temp, (1, n_samps))
        result[i,:] .= temp[1,:]
    end
    return reshape(result, size(input))
end

function apply_gate!(chip::Chip, tab::Vector{<:Integer}, bits::Vector{<:Integer})
    tensor_bmp = apply_noclean(chip.bitlines[bits])
    n_bits = length(bits)
    for (i,b) in zip(n_bits-1:-1:0, bits)
        bit_tab = tab .>> i .& 1
        chip.bitlines[b] = clean1(tensor_bmp, bit_tab)
    end
end

function minapply_gate!(chip::Chip, tab::Vector{<:Integer}, bits::Vector{<:Integer})
    sum_bmp, U = minapply(chip.bitlines[bits])
    n_bits = length(bits)
    for (i,b) in zip(n_bits-1:-1:0, bits)
        bit_tab = tab .>> i .& 1
        R = minapply_term(U, fill([0,1], length(bits)), bit_tab)
        chip.bitlines[b] = clean1_rl(sum_bmp, R)
    end
end

function apply_gate!(chip::Chip, g::ReversibleGate)
    apply_gate!(chip, g.perm, g.bits)
end

function apply_circuit!(chip::Chip, circuit::ReversibleCircuit)
    if circuit.n != length(chip.bitlines)
        throw(DimensionMismatch("Number of bitlines on the chip and the circuit don't match."))
    end
    for cg in circuit.gates
        apply_gate!(chip, cg)
    end
end

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

