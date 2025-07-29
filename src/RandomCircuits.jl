# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

using Random

function add_random_gates!(circuit::ReversibleCircuit, gate_bits)
    n_gates = length(gate_bits)
    sizehint!(circuit.gates, length(circuit.gates)+n_gates)
    for g_bits in gate_bits
        gsize = length(g_bits)
        perm = randperm(2^gsize) .- 1
        cg = ReversibleGate(perm, Vector{Int32}(g_bits))
        add_gate!(circuit, cg)
    end
end

function apply_random_gates!(chip::Chip, gate_bits)
    for g_bits in gate_bits
        gsize = length(g_bits)
        perm = randperm(2^gsize) .- 1
        apply_gate!(chip, perm, collect(g_bits))
    end
end

