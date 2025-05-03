function multiapply_noclean(
    chip::BareBMP,
    bits::Vector{<:Vector{<:Integer}}
)
    n = size(chip, 1)
    U = [(RSMInt(k), Vector{RSMInt}(bs)) for (k,bs) in enumerate(bits)]
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        A = Dict{Tuple{RSMInt, Vector{RSMInt}}, RSMInt}()
        S0 = [get!(A, (k, [chip[i,1].rows[j] for j in bs]), length(A)+1)
            for (k,bs) in U]
        S1 = [get!(A, (k, [chip[i,2].rows[j] for j in bs]), length(A)+1)
            for (k,bs) in U]
        mats[i,1] = RowSwitchMatrix(S0, length(A))
        mats[i,2] = RowSwitchMatrix(S1, length(A))
        U = Vector{Tuple{RSMInt, Vector{RSMInt}}}(undef, length(A))
        for k in keys(A)
            U[A[k]] = k
        end
    end
    return (mats, U)
end

"""
    multiapply(chip::BMP, bits::Vector{<:Vector{<:Integer}}, tabs::Vector{<:Vector{<:Integer}})

Uses a modified version of the direct-sum method to generate the joint BMP of
multiple functions in one sweep over the joint BMP `chip`. Each element of
`bits` is a vector indicating which of the output bits of `chip` enter into a
given function. The corresponding entry in `tabs` gives the truth table of that
function. For example, the call
```
multiapply(chip, [[1, 2], [1, 3]], [[0,0,0,1], [0,1,1,1]])
```
returns a BMP whose first output is the AND of the outputs `1` and `2` of
`chip`, and whose second output is the OR of the outputs `1` and `3`.
"""
function multiapply(
    chip::BMP,
    bits::Vector{<:Vector{<:Integer}},
    tabs::Vector{<:Vector{<:Integer}}
)
    mats, U = multiapply_noclean(chip.M, bits)
    R = fill(RSMInt(0), length(U))
    for (i, row) in enumerate(U)
        k, bs = row
        tabval = 0
        for bval in bs
            tabval = 2 * tabval + bval - 1
        end
        R[i] = tabs[k][tabval + 1]
    end
    return BMP(clean1_rl(mats, R), [0,1], copy(chip.order))
end

"""
    layerapply(bmp::BMP, bits::Vector{<:Vector{<:Integer}}, tabs::Vector{<:Vector{<:Integer}})

Given a joint BMP `chip` representing a reversible circuit, uses
[`multiapply`](@ref) to sythesize the BMP representing the application of
a layer of gates. Each element of `bits` is vector indicating the bits that
enter a gate, while the corresponding element of `tabs` indicates the
permutation implemented by the gate. For example, the call
```
layerapply(chip, [[1,4], [2,3]], [[0,1,3,2], [0,3,2,1]])
```
simulates the action of two CNOT gates: One on bitlines `1` and `4` where `1`
is the control, and another on bitlines `2` and `3` with `3` as the control.
"""
function layerapply(
    bmp::BMP,
    bits::Vector{<:Vector{<:Integer}},
    tabs::Vector{<:Vector{<:Integer}}
)
    n = length(bmp)
    bits_ = Vector{Vector{UInt8}}(undef, n)
    tabs_ = Vector{Vector{UInt8}}(undef, n)
    for (gbits, gtab) in zip(bits, tabs)
        gsize = length(gbits)
        for (s, bline) in zip(gsize-1:-1:0, gbits)
            bits_[bline] = gbits
            tabs_[bline] = (gtab .>> s) .& 1
        end
    end
    return multiapply(bmp, bits_, tabs_)
end
