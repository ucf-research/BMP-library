include("../src/TreeStructuredCircuit.jl")

using Random

function random_gates(q::Integer)
    all_gates = [shuffle(0:7) for i=1:3^q, j=1:q]
    return all_gates
end

function build_bmp(q::Integer, all_gates)
    n = 3^q
    bmp = BMP_join([BMP_bitline(i, n) for i in 1:n])
    for l=1:q
        bmp = TSC_apply_layer(q, l, bmp, all_gates[:,l])
    end
    return bmp
end

function build_circuit(q::Integer, all_gates)
    n = 3^q
    circuit = Circuit(n)
    for l=1:q
        TSC_apply_layer!(q, l, circuit, all_gates[:,l])
    end
    return circuit
end

function test(q::Integer, n_tests::Integer)
    n = 3^q
    all_gates = @time random_gates(q)
    circuit = @time build_circuit(q, all_gates)
    bmp = @time build_bmp(q, all_gates)
    tests_in = bitrand(n, n_tests)
    out1 = eval(circuit, tests_in)
    out2 = eval(bmp, tests_in)
    @show all(out1 .== out2)
end

let
    test(2, 100)
    test(2, 100)
    test(3, 1000)
    test(4, 1000)
    test(5, 1000)
    test(6, 1000)
    test(7, 1000)
end
