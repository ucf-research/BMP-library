# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

const RSMInt = Int32
struct RowSwitchMatrix
    rows::Vector{RSMInt}
    ncols::RSMInt
end

function RowSwitchMatrix(A::Matrix{<:Integer})
    N, M = size(A,1), size(A,2)
    cond = all((count(A.==1, dims=2) .== 1)
               .&& (count(A.==0, dims=2) .== M-1))
    if !cond
        throw(ArgumentError("Input argument is not row-switching."))
    end
    return RowSwitchMatrix([findfirst(A[i,:].==1) for i in 1:N], M)
end

function RowSwitchMatrix(n::Integer)
    return RowSwitchMatrix(Vector{RSMInt}(1:n), n)
end

function mult(a::RowSwitchMatrix, b::RowSwitchMatrix)
    if length(b.rows) != a.ncols
        throw(DimensionMismatch("Number of columns of a must equal the number of rows of b."))
    end
    return RowSwitchMatrix(b.rows[a.rows], b.ncols)
end

function mult_inplace(a::Vector{<:Integer}, b::RowSwitchMatrix)
    for i in eachindex(a)
        a[i] = b.rows[a[i]]
    end
end

function kron(mats::NTuple{N, RowSwitchMatrix}) where {N}
    nrows = prod(length(m.rows) for m in mats)
    cstrides = cumprod(m.ncols for m in Iterators.reverse(mats))
    S = Vector{RSMInt}(undef, nrows)
    data = reverse(map(m -> m.rows, mats))
    for (i, rs) in enumerate(Iterators.product(data...))
        col = rs[1] + sum(cstrides[i-1] * (rs[i]-1) for i=2:N)
        S[i] = col
    end
    return RowSwitchMatrix(S, cstrides[N])
end

function kron(mats::RowSwitchMatrix...)
    return kron(mats)
end

function kron(mats::Vector{RowSwitchMatrix})
    return kron(mats...)
end

function dsum(mats)
    nrows = sum(length(m.rows) for m in mats)
    ncols = sum(m.ncols for m in mats)
    array = fill(RSMInt(0), nrows)
    rstride = 0
    cstride = 0
    for m in mats
        mr = length(m.rows)
        for i=1:mr
            array[rstride+i] += cstride + m.rows[i]
        end
        rstride += mr
        cstride += m.ncols
    end
    return RowSwitchMatrix(array, ncols)
end

function bitmatrix(a::RowSwitchMatrix)
    nr = length(a.rows)
    nc = a.ncols
    result = BitMatrix(undef, (nr, nc))
    for i=1:nr, j=1:nc
        result[i,j] = a.rows[i] == j
    end
    return result
end

function to_matrix(a::RowSwitchMatrix)
    return [Int64(a.rows[i] == j) for i=1:length(a.rows), j=1:a.ncols]
end

function decompose(a::RowSwitchMatrix)
    u = RowSwitchMatrix(unique!(sort(a.rows)), a.ncols)
    s = RowSwitchMatrix(indexin(a.rows, u.rows), length(u.rows))
    return (s, u)
end
