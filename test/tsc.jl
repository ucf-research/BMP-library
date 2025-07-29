using BinaryMatrixProducts
using Test
using Random

@testset "Tree-structured circuit tests" begin
    for q=2:5
        chip = random_tsc_chip(q)
        @test volume(chip) < (3^q)^3
    end
    q = 3
    n_tests = 10
    for _ in 1:n_tests
        circuit = random_tsc(q)
        chip = Chip(circuit)
        test_inputs = bitrand(3^q, 1000)
        @test all(evalfunc(chip, test_inputs) .== evalfunc(circuit, test_inputs))
    end
end
