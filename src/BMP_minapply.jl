# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

function minapply_term(
    htab::AbstractArray,
    U::Dict{NTuple{N, RSMInt}, RSMInt},
    Rs::NTuple{N, <:AbstractArray}
) where {N}
    R = Vector{RSMInt}(undef, length(U))
    for (inds, i) in pairs(U)
        val = sum(2^(N-j) * Rs[j][inds[j]] for j=1:N)
        R[i] = htab[val+1]
    end
    return R
end

function minapply_noclean(bmps::NTuple{N, BareBMP}) where {N}
    n = size(bmps[1], 1)
    U = Dict{NTuple{N, RSMInt}, RSMInt}()
    U[ntuple(i -> RSMInt(1), N)] = 1
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        A = Dict{NTuple{N, RSMInt}, RSMInt}()
        sizehint!(A, 2 * length(U))
        S = [
            Vector{RSMInt}(undef, length(U)),
            Vector{RSMInt}(undef, length(U))
        ]
        for j=1:2
            for (rvals, rnum) in pairs(U)
                # rvals[x]: location of the non-zero column in the bmps[x] block on this row
                new_rvals = ntuple(x -> bmps[x][i,j].rows[rvals[x]], N)
                # new_rvals describes the new row obtained after matrix multiplication
                S[j][rnum] = get!(A, new_rvals, length(A)+1)
            end
        end
        for j=1:2
            mats[i,j] = RowSwitchMatrix(S[j], length(A))
        end
        U = A
    end
    return (mats, U)
end

function minapply_noclean(bmps::BareBMP...)
    return minapply_noclean(bmps)
end

function minapply_noclean(bmps::Array{BareBMP})
    return minapply_noclean(ntuple(i -> bmps[i], length(bmps)))
end

function minapply(
    htab::AbstractArray,
    bmps::NTuple{N, BareBMP},
    Rs::NTuple{N, <:AbstractArray}
) where {N}
    M, U = minapply_noclean(bmps)
    R = minapply_term(htab, U, Rs)
    return clean1_rl(M, R)
end

function minapply(
    htab::AbstractArray,
    bmps::Array{BareBMP},
    Rs::AbstractArray{<:AbstractArray}
)
    return minapply(htab, ntuple(i -> bmps[i], length(bmps)), ntuple(i -> Rs[i], length(Rs)))
end

function minapply_noclean(htab::AbstractArray, bmps::NTuple{N, BMP}) where {N}
    M, U = minapply_noclean(ntuple(i -> bmps[i].M, N))
    R = minapply_term(htab, U, ntuple(i -> bmps[i].R, N))
    return BMP(M, R, copy(bmps[1].order))
end

function minapply_noclean(htab::AbstractArray, bmps::BMP...)
    return minapply_noclean(htab, bmps)
end

function minapply_noclean(htab::AbstractArray, bmps::Array{BMP})
    return minapply_noclean(htab, ntuple(i -> bmps[i], length(bmps)))
end

"""
    minapply(htab::AbstractArray, bmps)

Implements the direct-sum APPLY operation. This creates a BMP for the Boolean
function ``h(f_1(\\vec{x}), \\dotsc, f_k(\\vec{x}))`` from the already known BMPs
of functions ``f_1,\\dotsc,f_n``.

The truth table of ``h(y_1,\\dotsc,y_k)`` must be given in `htab` such that
`htab[i]` the output of the function when its inputs ``y_1,\\dotsc,y_k``
are set to the bits of `i-1` going from most significant to the least, respectively.

The BMPs for ``f_1,\\dotsc,f_n`` must be given in `bmps`, either as a container
(such as a `NTuple` or `Vector`) or vararg of `BMP`. Where possible, prefer types
with sizes known at compile time.

See also [`apply`](@ref).
"""
function minapply(htab::AbstractArray, bmps::NTuple{N, BMP}) where {N}
    M, U = minapply_noclean(ntuple(i -> bmps[i].M, N))
    R = minapply_term(htab, U, ntuple(i -> bmps[i].R, N))
    return BMP(clean1_rl(M, R), RSMInt[0,1], copy(bmps[1].order))
end

function minapply(htab::AbstractArray, bmps::BMP...)
    return minapply(htab, bmps)
end

function minapply(htab::AbstractArray, bmps::Array{BMP})
    return minapply(htab, ntuple(i -> bmps[i], length(bmps)))
end
