const RSMInt = UInt32
struct RowSwitchMatrix
    rows::Vector{UInt32}
    ncols::UInt32
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

function RSM_mult(a::RowSwitchMatrix, b::RowSwitchMatrix)::RowSwitchMatrix
    if length(b.rows) != a.ncols
        throw(DimensionMismatch("Number of columns of a must equal the number of rows of b."))
    end
    return RowSwitchMatrix(b.rows[a.rows], b.ncols)
end

function RSM_kron(a::RowSwitchMatrix, b::RowSwitchMatrix)::RowSwitchMatrix
    Na = length(a.rows)
    Nb = length(b.rows)
    Ma = a.ncols
    Mb = b.ncols
    return RowSwitchMatrix(
        [(a.rows[k÷Nb+1]-1)*Mb + b.rows[k%Nb+1]
            for k in 0:(Na*Nb-1)],
        Ma * Mb
    )
end

function RSM_kron(mats::Vector{RowSwitchMatrix})::RowSwitchMatrix
    nrows = prod(length.(m.rows for m in mats))
    ncols = prod(m.ncols for m in mats)
    array = fill(RSMInt(0), nrows)
    stride = nrows
    for m in mats
        array .*= m.ncols
        mr = length(m.rows)
        stride = div(stride, mr) # distance between blocks
        cnt = div(nrows, stride) # number of blocks
        for i=0:cnt-1
            array[i*stride+1:i*stride+stride] .+= (m.rows[i % mr + 1] - 1)
        end
    end
    array .+= 1
    return RowSwitchMatrix(array, ncols)
end

function RSM_join(mats::Vector{RowSwitchMatrix})::RowSwitchMatrix
    nrows = sum(length.(m.rows for m in mats))
    ncols = sum(m.ncols for m in mats)
    array = fill(RSMInt(0), nrows)
    rstride = 0
    cstride = 0
    for m in mats
        mr = length(m.rows)
        array[rstride+1:rstride+mr] .+= cstride
        array[rstride+1:rstride+mr] .+= m.rows
        rstride += mr
        cstride += m.ncols
    end
    return RowSwitchMatrix(array, ncols)
end

function RSM_matrix(a::RowSwitchMatrix)
    return [Int64(a.rows[i] == j) for i=1:length(a.rows), j=1:a.ncols]
end

function SU_decomposition(a::RowSwitchMatrix)
    u = RowSwitchMatrix(unique!(sort(a.rows)), a.ncols)
    s = RowSwitchMatrix(indexin(a.rows, u.rows), length(u.rows))
    return (s, u)
end
