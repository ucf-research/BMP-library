include("../src/Chip.jl")
include("../src/Circuit.jl")

using Random

let
    n = 10
    n_gates = 10
    circuit = Circuit(n)
    chip = Chip(n)
    for gi=1:n_gates
        k = rand(2:4)
        bits = shuffle(1:n)
        perm = shuffle(0:2^k-1)
        add_gate!(circuit, perm, bits[1:k])
        Chip_apply!(chip, perm, bits[1:k])
    end
    bmp = Chip_join(chip)
    n_tests = 1000
    tests_in = bitrand(n, n_tests)
    tests1 = eval(circuit, tests_in)
    @show all(tests1 .== eval(bmp, tests_in))
    for it=1:20
        i = rand(1:n-1)
        BMP_swap!(bmp, i)
        @show i, Vector{Int64}(bmp.order)
        @show all(tests1 .== eval(bmp, tests_in))
    end
end

