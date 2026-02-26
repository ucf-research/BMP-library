using BinaryMatrixProducts
using BinaryMatrixProducts: max_dim
using Test
using Random

@testset "Reversible circuit operations" begin
    n = 12
    all_bits = collect(1:n)
    n_gates = 12
    n_tests = 10
    for _ in 1:n_tests
        circ = ReversibleCircuit(n)
        chip = Chip(n)
        for gi=1:n_gates
            shuffle!(all_bits)
            bits = all_bits[1:3]
            perm = randperm(2^length(bits)) .- 1
            add_gate!(circ, perm, bits)
            apply_gate!(chip, perm, bits)
        end
        test_inputs = bitrand(n, 1000)
        test_outputs = evalfunc(circ, test_inputs)
        # Circuit - Chip conversion correctness
        @test all(evalfunc(chip, test_inputs) .== test_outputs)
        chip1 = Chip(circ)
        @test all(evalfunc(chip1, test_inputs) .== test_outputs)
        # Circuit inversion correctness
        inv_circ = invert_circuit(circ)
        @test all(test_inputs .== evalfunc(inv_circ, test_outputs))
        apply_circuit!(chip1, inv_circ)
        @test check_equivalence(chip1, Chip(n))
    end
end
