using BinaryMatrixProducts
using Random
using Test

@testset "Gate conjugation" begin
    q = 2
    n = 3^q
    gate = ReversibleGate([0, 1, 3, 2], [div(n,2), div(n,2)+1])
    n_tests = 10
    for _ in 1:n_tests
        C = random_tsc(q)
        C1 = invert_circuit(C)
        gconj1 = Chip(n)
        minapply_circuit!(gconj1, C1)
        minapply_gate!(gconj1, gate)
        minapply_circuit!(gconj1, C)
        CC = conjugate_gate(gate, C)
        @test length(CC.gates) <= 2*length(C.gates)+1
        gconj2 = Chip(n)
        minapply_circuit!(gconj2, CC)
        @test check_equivalence(gconj1, gconj2)
    end
end
