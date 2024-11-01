include("BMP.jl")
include("Circuit.jl")

struct Chip
    bitlines::Vector{BareBMP}
    order::Vector{UInt32}
    position::Vector{UInt32}
    function Chip(n::Integer, order::Vector{<:Integer})
        position = fill(UInt32(0), length(order))
        position[order] .= order
        bitlines = [BMP_bitline_bare(xi, n) for xi in position]
        return new(bitlines, order, position)
    end
end

function Chip(n::Integer)
    return Chip(n, collect(1:n))
end

function eval(chip::Chip, input::BitArray)
    n = size(input, 1)
    n_samps = div(length(input), size(input, 1))
    result = BitArray(undef, (n, n_samps))
    for (i,bl) in enumerate(chip.bitlines)
        temp = eval(bl, input, [0,1], chip.order)
        temp = reshape(temp, (1, n_samps))
        result[i,:] .= temp[1,:]
    end
    return reshape(result, size(input))
end

function Chip_apply!(chip::Chip, tab::Vector{<:Integer}, bits::Vector{<:Integer})
    tensor_bmp = BMP_apply_noclean(chip.bitlines[bits])
    n_bits = length(bits)
    for (i,b) in zip(n_bits-1:-1:0, bits)
        bit_tab = tab .>> i .& 1
        chip.bitlines[b] = BMP_clean1(tensor_bmp, bit_tab)
    end
end

function Chip_minapply!(chip::Chip, tab::Vector{<:Integer}, bits::Vector{<:Integer})
    sum_bmp, U = BMP_minapply(chip.bitlines[bits])
    n_bits = length(bits)
    for (i,b) in zip(n_bits-1:-1:0, bits)
        bit_tab = tab .>> i .& 1
        R = BMP_minapply_R(U, fill([0,1], length(bits)), bit_tab)
        chip.bitlines[b] = BMP_clean1_rl(sum_bmp, R)
    end
end

function Chip_apply!(chip::Chip, g::CircuitGate)
    Chip_apply!(chip, g.tab, g.bits)
end

function Chip_apply!(chip::Chip, circuit::Circuit)
    if circuit.n != length(chip.bitlines)
        throw(DimensionMismatch("Number of bitlines on the chip and the circuit don't match."))
    end
    for cg in circuit.gates
        Chip_apply!(chip, cg)
    end
end

function Chip_join(chip::Chip)
    n = length(chip.bitlines)
    mats = Matrix{RowSwitchMatrix}(undef, (n, 2))
    for i=1:n, j=1:2
        mats[i,j] = RSM_join([bmp[i,j] for bmp in chip.bitlines])
    end
    R = fill(0, 2*n)
    R[2:2:end] .= 1
    return BMP_clean1(BMP(mats, R, copy(chip.order)))
end

function Chip_swap!(chip::Chip, i::Integer)
end

