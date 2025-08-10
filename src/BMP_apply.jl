# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

function apply_term(htab, Rs::NTuple{N, <:AbstractArray}) where {N}
    sz = prod(length.(Rs))
    R = Vector{RSMInt}(undef, sz)
    bvals = ntuple(i -> 2^(i-1), N)
    for (i, bits) in enumerate(Iterators.product(reverse(Rs)...))
        R[i] = htab[sum(bits[j] * bvals[j] for j=1:N) + 1]
    end
    return R
end

function apply_term(htab, Rs)
    return apply_term(htab, ntuple(i -> Rs[i], length(Rs)))
end

function apply_mats(bmps::NTuple{N, BareBMP}) where {N}
    n = size(bmps[1], 1)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        kron_mats = ntuple(x -> bmps[x][i,j], N)
        mats[i,j] = kron(kron_mats)
    end
    return mats
end

function apply_mats(bmps::BareBMP...)
    return apply_mats(bmps)
end

function apply_mats(bmps)
    return apply_mats(ntuple(i -> bmps[i], length(bmps)))
end

function apply(
    htab,
    bmps::NTuple{N, BareBMP},
    Rs::NTuple{N, <:AbstractArray};
    noclean::Bool=false
) where {N}
    M = apply_mats(bmps)
    R = apply_term(htab, Rs)
    if noclean
        return (M, R)
    end
    return (clean1(M, R), RSMInt[0,1])
end

function apply(
    htab,
    bmps::Array{BareBMP},
    Rs::Array{<:AbstractArray};
    noclean::Bool=false
)
    N = length(bmps)
    return apply(htab, ntuple(i -> bmps[i], N), ntuple(i -> Rs[i], N); noclean=noclean)
end

"""
    apply(htab, bmps; noclean::Bool=false)

Implements the direct-product APPLY operation. This creates a BMP for the Boolean
function ``h(f_1(\\vec{x}), \\dotsc, f_k(\\vec{x}))`` from the already known BMPs
of functions ``f_1,\\dotsc,f_n``.

The truth table of ``h(y_1,\\dotsc,y_k)`` must be given in `htab` such that
`htab[i]` the output of the function when its inputs ``y_1,\\dotsc,y_k``
are set to the bits of `i-1` going from most significant to the least, respectively.

The BMPs for ``f_1,\\dotsc,f_n`` must be given in `bmps`, either as a container
(such as a `NTuple` or `Vector`) or vararg of `BMP`. Where possible, prefer types
with sizes known at compile time.

By default, the CLEAN operation is invoked to reduce the resulting BMP to its
canonical form. If this is not desired, the keyword argument `noclean` should be
set to `true`.

See also [`minapply`](@ref).
"""
function apply(htab, bmps::NTuple{N, BMP}; noclean::Bool=false) where {N}
    mats = ntuple(i -> bmps[i].M, N)
    Rs = ntuple(i -> bmps[i].R, N)
    M = apply_mats(mats)
    R = apply_term(htab, Rs)
    if noclean
        return BMP(M, R, copy(bmps[1].order))
    end
    return BMP(clean1(M, R), RSMInt[0,1], copy(bmps[1].order))
end

function apply(htab, bmps::BMP...; noclean::Bool=false)
    return apply(htab, bmps; noclean=noclean)
end

function apply(htab, bmps::Array{BMP}; noclean::Bool=false)
    N = length(bmps)
    return apply(htab, ntuple(i -> bmps[i], N); noclean=noclean)
end
