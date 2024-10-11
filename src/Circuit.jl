struct CircuitGate
    tab::Vector{UInt16}
    bits::Vector{UInt32}
    function CircuitGate(tab::Vector{<:Integer}, bits::Vector{<:Integer})
        # Enforce the number of bits
        if length(tab) != 2^length(bits)
            throw(DimensionMismatch("The truth table of the gate does not have the correct number of entries."))
        end
        return new(tab, bits)
    end
end

struct Circuit
    n::UInt32
    gates::Vector{CircuitGate}
    function Circuit(n::Integer)
        return new(n, CircuitGate[])
    end
end

function eval(circuit::Circuit, input::BitArray)
    n = size(input, 1)
    if n != circuit.n
        throw(DimensionMismatch("The circuit input length and the number of circuit registers don't match."))
    end
    input_ = reshape(copy(input), (n, div(length(input), n)))
    for j=1:size(input_, 2)
        for g in circuit.gates
            val = 0
            for b in g.bits
                val = 2 * val + input_[b,j]
            end
            output = g.tab[val+1]
            for (k,b) in zip(length(g.bits)-1:-1:0, g.bits)
                input_[b,j] = output >> k & 1
            end
        end
    end
    return reshape(input_, size(input))
end

function add_gate!(circuit::Circuit, tab::Vector{<:Integer}, bits::Vector{<:Integer})
    push!(circuit.gates, CircuitGate(tab, bits))
end

