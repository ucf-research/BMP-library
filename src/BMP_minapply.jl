# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

function minapply_term(
    U::Vector{Tuple{RSMInt, RSMInt}},
    R1::Vector{<:Integer},
    R2::Vector{<:Integer},
    htab::Vector{<:Integer}
)
    R = fill(RSMInt(0), length(U))
    for (i,pair) in enumerate(U)
        i1, i2 = pair
        val = 2 * R1[i1] + R2[i2]
        R[i] = htab[val+1]
    end
    return R
end

function minapply_noclean(
    bmp1::BareBMP,
    bmp2::BareBMP
)
    n = size(bmp1,1)
    U = [(RSMInt(1),RSMInt(1))]
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i in axes(bmp1,1)
        A = Dict{Tuple{RSMInt, RSMInt}, RSMInt}()
        sizehint!(A, 2 * length(U))
        S0 = [get!(A, (bmp1[i,1].rows[i1], bmp2[i,1].rows[i2]), length(A)+1)
            for (i1, i2) in U]
        S1 = [get!(A, (bmp1[i,2].rows[i1], bmp2[i,2].rows[i2]), length(A)+1)
            for (i1, i2) in U]
        mats[i,1] = RowSwitchMatrix(S0, length(A))
        mats[i,2] = RowSwitchMatrix(S1, length(A))
        U = fill((RSMInt(0), RSMInt(0)), length(A))
        for k in keys(A)
            U[A[k]] = k
        end
    end
    return (mats, U)
end

function minapply(
    bmp1::BareBMP,
    bmp2::BareBMP,
    R1::Vector{<:Integer},
    R2::Vector{<:Integer},
    htab::Vector{<:Integer}
)
    M, U = minapply_noclean(bmp1, bmp2)
    R = minapply_term(U, R1, R2, htab)
    return clean1_rl(M, R)
end

function minapply_noclean(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
    M, U = minapply_noclean(bmp1.M, bmp2.M)
    return BMP(M, minapply_term(U, bmp1.R, bmp2.R, htab), copy(bmp1.order))
end

"""
    minapply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})

Returns the result of the direct-sum APPLY method performed on `bmp1` and
`bmp2`. The truth table of the function is given in `htab` in the order where
`bmp1` value is the most significant bit.

See also [`apply`](@ref).
"""
function minapply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})::BMP
    M, U = minapply_noclean(bmp1.M, bmp2.M)
    R = minapply_term(U, bmp1.R, bmp2.R, htab)
    return BMP(clean1_rl(M, R), RSMInt[0,1], copy(bmp1.order))
end

function minapply_noclean(bmps::Vector{BareBMP})
    n = size(bmps[1], 1)
    k = length(bmps)
    init_view = CustomView(fill(RSMInt(1), k), 1, k)
    U = Dict{CustomView, RSMInt}()
    U[init_view] = 1
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    site_mats = Matrix{RowSwitchMatrix}(undef, (k,2))
    for i=1:n
        for b=1:2, j=1:k
            site_mats[j,b] = bmps[j][i,b]
        end
        chi = length(U)
        A = propagate_mat(U, site_mats)
        #
        U = Dict{CustomView, RSMInt}()
        sizehint!(U, 2 * chi)
        S = [fill(RSMInt(0), chi), fill(RSMInt(0), chi)]
        for b=0:1
            for r=1:chi
                rbegin = A.rbegin[r + b * chi]
                rend = A.rend[r + b * chi]
                row_key = CustomView(A.data, rbegin, rend)
                S[b+1][r] = get!(U, row_key, length(U)+1)
            end
        end
        mats[i,1] = RowSwitchMatrix(S[1], length(U))
        mats[i,2] = RowSwitchMatrix(S[2], length(U))
    end
    return (mats, U)
end

function minapply_term(
    U::Dict{CustomView, RSMInt},
    Rs::Vector{<:Vector{<:Integer}},
    htab::Vector{<:Integer}
)
    R = fill(RSMInt(0), length(U))
    for (vw, k) in pairs(U)
        val = 0
        len = vw.re - vw.rb + 1
        for (i,j) in zip(vw.data[vw.rb:vw.re], 1:len)
            val = 2 * val + Rs[j][i]
        end
        R[k] = htab[val+1]
    end
    return R
end

function minapply(
    bmps::Vector{BareBMP},
    Rs::Vector{<:Vector{<:Integer}},
    htab::Vector{<:Integer}
)
    mats, U = minapply_noclean(bmps)
    R = minapply_term(U, Rs, htab)
    return clean1_rl(mats, R)
end

function minapply_noclean(bmps::Vector{BMP}, htab::Vector{<:Integer})
    mats, U = minapply_noclean([bmp.M for bmp in bmps])
    R = minapply_term(U, [bmp.R for bmp in bmps], htab)
    return BMP(mats, R, copy(bmps[1].order))
end

"""
    minapply(bmps::Vector{BMP}, htab::Vector{<:Integer})

Returns the result of the direct-sum APPLY method performed on the elements
of `bmps`. The truth table of the function is given in `htab`, in the order
where the most significant bit corresponds to `bmps[1]`.

See also [`apply`](@ref).
"""
function minapply(bmps::Vector{BMP}, htab::Vector{<:Integer})
    mats, U = minapply_noclean([bmp.M for bmp in bmps])
    R = minapply_term(U, [bmp.R for bmp in bmps], htab)
    mats_ = clean1_rl(mats, R)
    return BMP(mats_, RSMInt[0,1], copy(bmps[1].order))
end

