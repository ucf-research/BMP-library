# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

function apply_term(R1::Vector{<:Integer}, R2::Vector{<:Integer}, htab::Vector{<:Integer})
    R = [htab[2*i+j+1] for (j,i) in Iterators.product(R2, R1)]
    R = reshape(R, length(R))
    return R
end

function apply_term(Rs::Vector{<:Vector{<:Integer}}, htab::Vector{<:Integer})
    nrows = prod(length.(R for R in Rs))
    inds = fill(0, nrows)
    stride = nrows
    for R in Rs
        inds .*= 2
        mr = length(R)
        stride = div(stride, mr)
        cnt = div(nrows, stride)
        for i=0:cnt-1
            inds[i*stride+1:i*stride+stride] .+= R[i % mr + 1]
        end
    end
    return htab[inds .+ 1]
end

function apply_noclean(bmp1::BareBMP, bmp2::BareBMP)::BareBMP
    n = size(bmp1, 1)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        mats[i,j] = kron(bmp1[i,j], bmp2[i,j])
    end
    return mats
end

function apply(bmp1::BareBMP, bmp2::BareBMP, htab::Vector{<:Integer})
    # Assume canonical R
    return clean1(apply_noclean(bmp1, bmp2), htab)
end

function apply(bmp1::BareBMP, bmp2::BareBMP, R1::Vector{<:Integer}, R2::Vector{<:Integer}, htab::Vector{<:Integer})
    return clean1(apply_noclean(bmp1, bmp2), apply_term(R1, R2, htab))
end

function apply_noclean(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
    mats = apply_noclean(bmp1.M, bmp2.M)
    return BMP(mats, apply_term(bmp1.R, bmp2.R, htab), copy(bmp1.order))
end

"""
    apply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})

Returns the result of the direct-product APPLY method performed on `bmp1` and
`bmp2`. The truth table of the function is given in `htab` in the order where
`bmp1` value is the most significant bit.

See also [`minapply`](@ref).
"""
function apply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
    return clean1(apply_noclean(bmp1, bmp2, htab))
end

function apply_noclean(bmps::Vector{BareBMP})::BareBMP
    n = size(bmps[1], 1)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        mats[i,j] = kron([M[i,j] for M in bmps])
    end
    return mats
end

function apply(bmps::Vector{BareBMP}, htab::Vector{<:Integer})
    return clean1(apply_noclean(bmps), htab)
end

function apply(bmps::Vector{BareBMP}, Rs::Vector{<:Vector{<:Integer}}, htab::Vector{<:Integer})
    return clean1(apply_noclean(bmps), apply_term(Rs, htab))
end

function apply_noclean(bmps::Vector{BMP}, htab::Vector{<:Integer})
    mats = apply_noclean([bmp.M for bmp in bmps])
    return BMP(mats, apply_term([bmp.R for bmp in bmps], htab), copy(bmps[1].order))
end

"""
    apply(bmps::Vector{BMP}, htab::Vector{<:Integer})

Returns the result of the direct-product APPLY method performed on the elements
of `bmps`. The truth table of the function is given in `htab`, in the order
where the most significant bit corresponds to `bmps[1]`.

See also [`minapply`](@ref)
"""
function apply(bmps::Vector{BMP}, htab::Vector{<:Integer})
    return clean1(apply_noclean(bmps, htab))
end
