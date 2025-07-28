# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

struct ReversibleGate
    perm::Vector{UInt32}
    bits::Vector{Int32}
    function ReversibleGate(perm::Vector{<:Integer}, bits::Vector{<:Integer})
        # Enforce the number of bits
        if length(perm) != 2^length(bits)
            throw(DimensionMismatch("The permutation table of the gate does not have the correct number of entries."))
        end
        return new(perm, bits)
    end
end

function invert_gate(gate::ReversibleGate)
    iperm = copy(gate.perm)
    for (i, val) in enumerate(gate.perm)
        iperm[val+1] = i-1
    end
    ibits = copy(gate.bits)
    return ReversibleGate(iperm, ibits)
end

struct ReversibleCircuit
    n::Int32
    gates::Vector{ReversibleGate}
    function ReversibleCircuit(n::Integer)
        return new(n, ReversibleGate[])
    end
    function ReversibleCircuit(n::Integer, gates::Vector{ReversibleGate})
        return new(n, gates)
    end
end

function evalfunc(circuit::ReversibleCircuit, input::AbstractArray)
    n = size(input, 1)
    if n != circuit.n
        throw(DimensionMismatch("Input length does not match circuit register count."))
    end
    input_ = reshape(copy(input), (n, :))
    for j in axes(input_, 2)
        for g in circuit.gates
            val = 0
            for b in g.bits
                val = 2 * val + input_[b,j]
            end
            output = g.perm[val+1]
            for (k,b) in zip(length(g.bits)-1:-1:0, g.bits)
                input_[b,j] = output >> k & 1
            end
        end
    end
    return reshape(input_, size(input))
end

function add_gate!(circuit::ReversibleCircuit, perm::Vector{<:Integer}, bits::Vector{<:Integer})
    push!(circuit.gates, ReversibleGate(perm, bits))
end

function add_gate!(circuit::ReversibleCircuit, gate::ReversibleGate)
    push!(circuit.gates, gate)
end

function invert_circuit(circuit::ReversibleCircuit)
    n_gates = length(circuit.gates)
    gates = similar(circuit.gates, n_gates)
    for (i,g) in enumerate(Iterators.reverse(circuit.gates))
        gates[i] = invert_gate(g)
    end
    return ReversibleCircuit(circuit.n, gates)
end
