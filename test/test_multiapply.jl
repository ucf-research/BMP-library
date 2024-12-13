include("../src/BMP.jl")
include("../src/Circuit.jl")

using Random

function build_bmp(bmp, all_bits, all_gates)
    n_layers = length(all_bits)
    bmp_ = bmp
    for i in 1:n_layers
        bmp_ = @time BMP_layerapply(bmp_, all_bits[i], all_gates[i])
    end
    return bmp_
end

function test_generation(q, k, n_layers, n_tests)
    n = k * q
    circuit = Circuit(n)
    println("q = $q, k = $k, # layers = $n_layers, # tests = $n_tests")
    all_bits = []
    all_gates = []
    for _ in 1:n_layers
        bitperm = shuffle(1:n)
        bits = Vector{UInt8}[]
        tabs = Vector{UInt8}[]
        for i=1:q
            b = bitperm[k*(i-1)+1:k*i]
            t = shuffle(0:2^k-1)
            add_gate!(circuit, t, b)
            push!(bits, b)
            push!(tabs, t)
        end
        push!(all_bits, bits)
        push!(all_gates, tabs)
    end
    bmp_ = BMP_join([BMP_bitline(i, n) for i in 1:n])
    bmp = @time build_bmp(bmp_, all_bits, all_gates)
    inputs = bitrand(n, n_tests)
    outs1 = eval(circuit, inputs)
    outs2 = eval(bmp, inputs)
    @show all(outs1 .== outs2)
end

let
    test_generation(3, 3, 5, 100)
    test_generation(3, 3, 5, 100)
    test_generation(4, 3, 5, 100)
    test_generation(5, 3, 5, 100)
    test_generation(6, 3, 6, 100)
    test_generation(7, 3, 7, 100)
    test_generation(8, 3, 8, 100)
end
