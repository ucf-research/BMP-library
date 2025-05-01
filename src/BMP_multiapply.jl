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
