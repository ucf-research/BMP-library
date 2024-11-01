include("../src/Chip.jl")

function build_chip0(n, all_bits, all_gates)
    chip = Chip(n)
    for (bits, gate) in zip(all_bits, all_gates)
        Chip_apply!(chip, gate, bits)
    end
    return chip
end

function build_chip1(n, all_bits, all_gates)
    chip = Chip(n)
    for (bits, gate) in zip(all_bits, all_gates)
        Chip_minapply!(chip, gate, bits)
    end
    return chip
end

using Random

let
    # Run the functions first so that they get compiled
    build_chip0(3, [[1,2,3]], [[0,1,2,3,4,5,6,7]])
    build_chip1(3, [[1,2,3]], [[0,1,2,3,4,5,6,7]])
    # Randomly generate a number of gates
    n = 10
    n_gates = 2*n
    all_bits = []
    all_gates = []
    for i=1:n_gates
        k = rand(2:4)
        bits = shuffle(1:n)[1:k]
        gate = shuffle(0:2^k-1)
        push!(all_bits, bits)
        push!(all_gates, gate)
    end
    # Build two chips using the two different algorithms for APPLY and compare performance
    chip0 = @time build_chip0(n, all_bits, all_gates)
    chip1 = @time build_chip1(n, all_bits, all_gates)
    # Ensure that the two chips produce the same outputs
    bmp0 = Chip_join(chip0)
    bmp1 = Chip_join(chip1)
    tests = bitrand(n, 1000)
    @show all(eval(bmp0, tests) .== eval(bmp1, tests))
end
