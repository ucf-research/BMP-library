# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

function apply_term(htab::AbstractArray, Rs::NTuple{N, <:AbstractArray}) where {N}
    sz = prod(length.(Rs))
    R = Vector{RSMInt}(undef, sz)
    bvals = ntuple(i -> 2^(i-1), N)
    for (i, bits) in enumerate(Iterators.product(reverse(Rs)...))
        R[i] = htab[sum(bits[j] * bvals[j] for j=1:N) + 1]
    end
    return R
end

function apply_noclean(bmps::NTuple{N, BareBMP}) where {N}
    n = size(bmps[1], 1)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        kron_mats = ntuple(x -> bmps[x][i,j], N)
        mats[i,j] = kron(kron_mats)
    end
    return mats
end

function apply_noclean(bmps::BareBMP...)
    apply_noclean(bmps)
end

function apply_noclean(bmps::Array{BareBMP})
    N = length(bmps)
    return apply_noclean(ntuple(i -> bmps[i], N))
end

function apply(htab::AbstractArray, bmps::NTuple{N, BareBMP}, Rs::NTuple{N, <:AbstractArray}) where {N}
    return clean1(apply_noclean(bmps), apply_term(htab, Rs))
end

function apply(htab::AbstractArray, bmps::Array{BareBMP}, Rs::AbstractArray{<:AbstractArray})
    N = length(bmps)
    return apply(htab, ntuple(i -> bmps[i], N), ntuple(i -> Rs[i], N))
end

function apply_noclean(htab::AbstractArray, bmps::NTuple{N, BMP}) where {N}
    mats = ntuple(i -> bmps[i].M, N)
    Rs = ntuple(i -> bmps[i].R, N)
    return BMP(apply_noclean(mats), apply_term(htab, Rs), copy(bmps[1].order))
end

function apply_noclean(htab::AbstractArray, bmps::BMP...)
    return apply_noclean(htab, bmps)
end

function apply_noclean(htab::AbstractArray, bmps::Array{BMP})
    N = length(bmps)
    return apply_noclean(htab, ntuple(i -> bmps[i], N))
end

"""
    apply(htab::AbstractArray, bmps)

Implements the direct-product APPLY operation. This creates a BMP for the Boolean
function ``h(f_1(\\vec{x}), \\dotsc, f_k(\\vec{x}))`` from the already known BMPs
of functions ``f_1,\\dotsc,f_n``.

The truth table of ``h(y_1,\\dotsc,y_k)`` must be given in `htab` such that
`htab[i]` the output of the function when its inputs ``y_1,\\dotsc,y_k``
are set to the bits of `i-1` going from most significant to the least, respectively.

The BMPs for ``f_1,\\dotsc,f_n`` must be given in `bmps`, either as a container
(such as a `NTuple` or `Vector`) or vararg of `BMP`. Where possible, prefer types
with sizes known at compile time.

See also [`minapply`](@ref).
"""
function apply(htab::AbstractArray, bmps::NTuple{N, BMP}) where {N}
    mats = ntuple(i -> bmps[i].M, N)
    Rs = ntuple(i -> bmps[i].R, N)
    return BMP(
        clean1(apply_noclean(mats), apply_term(htab, Rs)),
        RSMInt[0,1],
        copy(bmps[1].order)
    )
end

function apply(htab::AbstractArray, bmps::BMP...)
    return apply(htab, bmps)
end

function apply(htab::AbstractArray, bmps::Array{BMP})
    N = length(bmps)
    return apply(htab, ntuple(i -> bmps[i], N))
end
