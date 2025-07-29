# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

function tsc_ordering(q::Integer, layers::Vector{<:Integer}; k::Integer=3)
    result = fill(0, k^q)
    terms = k .^(layers .- 1)
    for bit=0:k^q-1
        z = [div(bit, k^l) % k for l=0:q-1]
        loc = sum(terms .* z)
        result[loc+1] = bit
    end
    return result
end

struct TSCLayerGates
    n_gates::Int32
    n_clgates::Int32
    gsize::Int32
    function TSCLayerGates(q::Integer, l::Integer; k::Integer=3)
        l = 1 + (l-1) % q
        n_gates = k^(q-1)
        n_clgates = k^(l-1)
        new(n_gates, n_clgates, k)
    end
end

function Base.length(it::TSCLayerGates)
    return it.n_gates
end

function Base.iterate(it::TSCLayerGates, state=(1,0,0))
    gi, clstart, i = state
    if gi > it.n_gates
        return nothing
    end
    fbit = clstart + i + 1
    lbit = clstart + i + 1 + (it.gsize - 1) * it.n_clgates
    return (
        fbit:it.n_clgates:lbit,
        (
            gi+1,
            i == it.n_clgates-1 ? clstart + it.gsize*it.n_clgates : clstart,
            i == it.n_clgates-1 ? 0 : i+1
        )
    )
end

struct TSCGates
    n_gates::Int32
    n_layers::Int32
    gsize::Int32
    layers::Vector{TSCLayerGates}
    function TSCGates(q::Integer, layers; k::Integer=3)
        n_layers = length(layers)
        n_gates = n_layers * 3^(q-1)
        layer_iters = [TSCLayerGates(q, l; k) for l in layers]
        return new(n_gates, n_layers, k, layer_iters)
    end
end

function Base.length(it::TSCGates)
    return it.n_gates
end

function Base.iterate(it::TSCGates, state=(1, (1,0,0)))
    i, substate = state
    iter_output = iterate(it.layers[i], substate)
    while isnothing(iter_output)
        if i == it.n_layers
            return nothing
        end
        i = i+1
        iter_output = iterate(it.layers[i])
    end
    output, next_substate = iter_output
    return (output, (i, next_substate))
end

function tsc_gate_inputs(q::Integer, layers; k::Integer=3)
    return TSCGates(q, layers; k)
end

function tsc_gate_inputs(q::Integer, lmin::Integer, lmax::Integer; k::Integer=3)
    return tsc_gate_inputs(q, lmin:lmax; k)
end

function tsc_gate_inputs(q::Integer, lmax::Integer; k::Integer=3)
    return tsc_gate_inputs(q, 1:lmax; k)
end

function tsc_gate_inputs(q::Integer; k::Integer=3)
    return tsc_gate_inputs(q, 1:q; k)
end

function random_tsc(q::Integer, layers; k::Integer=3)
    circuit = ReversibleCircuit(3^q)
    add_random_gates!(circuit, tsc_gate_inputs(q, layers; k))
    return circuit
end

function random_tsc(q::Integer, lmin::Integer, lmax::Integer; k::Integer=3)
    return random_tsc(q, lmin:lmax; k)
end

function random_tsc(q::Integer, lmax::Integer; k::Integer=3)
    return random_tsc(q, 1:lmax; k)
end

function random_tsc(q::Integer; k::Integer=3)
    return random_tsc(q, 1:q; k)
end

function random_tsc_chip(q::Integer, layers; k::Integer=3)
    chip = Chip(3^q)
    apply_random_gates!(chip, tsc_gate_inputs(q, layers; k))
    return chip
end

function random_tsc_chip(q::Integer, lmin::Integer, lmax::Integer; k::Integer=3)
    return random_tsc_chip(q, lmin:lmax; k)
end

function random_tsc_chip(q::Integer, lmax::Integer; k::Integer=3)
    return random_tsc_chip(q, 1:lmax; k)
end

function random_tsc_chip(q::Integer; k::Integer=3)
    return random_tsc_chip(q, 1:q; k)
end
