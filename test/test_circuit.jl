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
        @time add_gate!(circuit, perm, bits[1:k])
        @time Chip_apply!(chip, perm, bits[1:k])
    end
    bmp = Chip_join(chip)
    n_tests = 1000
    tests_in = bitrand(n, n_tests)
    tests1 = eval(circuit, tests_in)
    tests2 = eval(chip, tests_in)
    @show all(tests1 .== tests2)
end
