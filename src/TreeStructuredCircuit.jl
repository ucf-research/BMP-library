# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

function effective_order(q::Integer, layers::Vector{<:Integer}, k::Integer=3)
    result = fill(0, k^q)
    terms = k .^(layers .- 1)
    for bit=0:k^q-1
        z = [div(bit, k^l) % k for l=0:q-1]
        loc = sum(terms .* z)
        result[loc+1] = bit
    end
    return result
end

struct TSCGates
    n_gates::UInt32
    n_clgates::UInt32
    gsize::UInt32
    function TSCGates(q::Integer, l::Integer, k::Integer=3)
        l = 1 + (l-1) % q
        n_gates = k^(q-1)
        n_clgates = k^(l-1)
        new(n_gates, n_clgates, k)
    end
end

function Base.iterate(it::TSCGates, state=(1,0,0))
    gi, clstart, i = state
    if gi > it.n_gates
        return nothing
    end
    fbit = clstart + i + 1
    lbit = clstart + i + 1 + (it.gsize - 1) * it.n_clgates
    return (
        (gi, fbit:it.n_clgates:lbit),
        (
            gi+1,
            i == it.n_clgates-1 ? clstart + it.gsize*it.n_clgates : clstart,
            i == it.n_clgates-1 ? 0 : i+1
        )
    )
end

function TSC_apply_layer!(
    q::Integer,
    l::Integer,
    circuit::Circuit,
    gates::Vector{<:Vector{<:Integer}},
    k::Integer=3
)
    n = circuit.n
    if n != k^q
        throw(ArgumentError("The BMP does not have the required number of reqisters."))
    end
    for (gi, bs) in TSCGates(q, l, k)
        gbits = collect(bs)
        add_gate!(circuit, gates[gi], gbits)
    end
end

function TSC_apply_layer(
    q::Integer,
    l::Integer,
    bmp::BMP,
    gates::Vector{<:Vector{<:Integer}},
    k::Integer=3
)
    n = length(bmp)
    if n != k^q
        throw(ArgumentError("The BMP does not have the required number of reqisters."))
    end
    bits_ = Vector{Vector{UInt16}}(undef, n)
    tabs_ = Vector{Vector{UInt16}}(undef, n)
    for (gi, bs) in TSCGates(q, l, k)
        gbits = collect(bs)
        for (s, bline) in zip(k-1:-1:0, gbits)
            bits_[bline] = gbits
            tabs_[bline] = (gates[gi] .>> s) .& 1
        end
    end
    return BMP_multiapply(bmp, bits_, tabs_)
end
