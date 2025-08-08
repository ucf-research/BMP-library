# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

"""
    ReversibleGate

A data type that encapsulates a gate in a reversible circuit. The `bits` field indicates
the bitlines that the gate acts on, while `perm` contains the permutation implemented
on those bitlines.
"""
struct ReversibleGate
    perm::Vector{UInt32}
    bits::Vector{Int32}
    function ReversibleGate(perm, bits)
        # Enforce the number of bits
        if length(perm) != 2^length(bits)
            throw(DimensionMismatch("The permutation table of the gate does not have the correct number of entries."))
        end
        return new(perm, bits)
    end
end

"""
    invert_gate(gate::ReversibleGate)

Returns the inverse of `gate`. This is a gate that acts on the same bitlines as `gate`,
but with the inverse of the permutation implemented by `gate`.
"""
function invert_gate(gate::ReversibleGate)
    iperm = copy(gate.perm)
    for (i, val) in enumerate(gate.perm)
        iperm[val+1] = i-1
    end
    ibits = copy(gate.bits)
    return ReversibleGate(iperm, ibits)
end

"""
    ReversibleGate

A data type that encapsulates a reversible circuit. This is essentially an ordered
collection of `ReversibleGate`s.
"""
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

"""
    evalfunc(circuit::ReversibleCircuit, input::AbstractArray)

Evaluates the circuit for the inputs given in `input`. This is done gate by gate.
`input` specifies the input vectors following the same convention as in `evalfunc` for
`BMP`, such that a slice (:,i1,i2,...,ik) is one input set.
"""
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

"""
    add_gate!(circuit::ReversibleGate, perm, bits)

Adds a gate specified by a permutation `perm` and bitlines `bits` to the gate. The added
gate acts after all the existing gates in the circuit.
"""
function add_gate!(circuit::ReversibleCircuit, perm, bits)
    push!(circuit.gates, ReversibleGate(perm, bits))
end

"""
    add_gate!(circuit::ReversibleCircuit, gate::ReversibleGate)

Adds a gate specified by a permutation `perm` and bitlines `bits` to the gate. The added
gate acts after all the existing gates in the circuit.
"""
function add_gate!(circuit::ReversibleCircuit, gate::ReversibleGate)
    push!(circuit.gates, gate)
end

"""
    invert_circuit(circuit::ReversibleCircuit)

Returns the inversion of `circuit`.
"""
function invert_circuit(circuit::ReversibleCircuit)
    n_gates = length(circuit.gates)
    gates = similar(circuit.gates, n_gates)
    for (i,g) in enumerate(Iterators.reverse(circuit.gates))
        gates[i] = invert_gate(g)
    end
    return ReversibleCircuit(circuit.n, gates)
end
