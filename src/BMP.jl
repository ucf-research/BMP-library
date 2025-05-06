# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

"""
    BareBMP

An alias for `Matrix{RowSwitchMatrix}`. Describes matrix product without
the terminal vector or the variable ordering information, for cases where
these should be handled separately.
"""
const BareBMP = Matrix{RowSwitchMatrix}

"""
    BMP

The main type used to represent binary matrix products. In addition to the
matrices, contains fields for the terminal vector and variable ordering.
"""
struct BMP
    M::Matrix{RowSwitchMatrix}
    R::Vector{RSMInt}
    order::Vector{UInt32}
    position::Vector{UInt32}
    function BMP(M::Matrix{RowSwitchMatrix}, R::Vector{<:Integer}, order::Vector{<:Integer})
        position = fill(UInt32(0), length(order))
        # position is the inverse of the permutation given in order
        position[order] .= order
        return new(M, R, order, position)
    end
end

# Constructors and initializers
function bare_bmp(val::Integer, n::Integer)::BareBMP
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n-1
        mats[i,:] = [RowSwitchMatrix(1), RowSwitchMatrix(1)]
    end
    mats[n,1] = RowSwitchMatrix(RSMInt[val+1], 2)
    mats[n,2] = RowSwitchMatrix(RSMInt[val+1], 2)
    return mats
end

"""
    BMP(val::Integer, n::Integer)

Returns a BMP for the constant function of value `val` of `n` input bits.
"""
function BMP(val::Integer, n::Integer)
    mats = bare_bmp(val, n)
    return BMP(mats, [0,1], collect(1:n))
end

"""
    BMP(val::Integer, order::Vector{<:Integer})

Returns a BMP for the constant function of value `val` with the variable
ordering `order`.
"""
function BMP(val::Integer, order::Vector{<:Integer})
    mats = bare_bmp(val, length(order))
    return BMP(mats, [0,1], copy(order))
end

function projbmp_bare(xi::Integer, n::Integer)::BareBMP
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:xi-1
        mats[i,:] = [RowSwitchMatrix(1), RowSwitchMatrix(1)]
        i += 1
    end
    mats[xi,1] = RowSwitchMatrix(RSMInt[1], 2)
    mats[xi,2] = RowSwitchMatrix(RSMInt[2], 2)
    for i=xi+1:n
        mats[i,:] = [RowSwitchMatrix(2), RowSwitchMatrix(2)]
    end
    return mats
end

"""
    projbmp(xi::Integer, n::Integer)

Returns a BMP of the projection to `xi` of `n` input variables,
i.e. ``f(x_1,\\dotsc, x_i, \\dotsc, x_n) = x_i``.
"""
function projbmp(xi::Integer, n::Integer)
    mats = projbmp_bare(xi, n)
    return BMP(mats, [0,1], collect(1:n))
end

"""
    projbmp(xi::Integer, order::Vector{<:Integer})

Returns the BMP of the projection to `xi` for the variable ordering `order`.
"""
function projbmp(xi::Integer, order::Vector{<:Integer})
    pos = findfirst(order .== xi)
    mats = projbmp_bare(pos, length(order))
    return BMP(mats, [0,1], copy(order))
end

"""
    Base.length(bmp::BMP)

Returns the number of input bits of `bmp`.
"""
function Base.length(bmp::BMP)
    return size(bmp.M, 1)
end

function bonddims(bmp::BareBMP)
    return [length(m.rows) for m in bmp[:,1]]
end

function bonddims(bmp::BareBMP, i::Integer)
    return length(bmp[i,1].rows)
end

"""
    bonddims(bmp::BMP)

Returns the bond dimensions of the BMP as an array. The terminal vector size
is not included in this array.
"""
function bonddims(bmp::BMP)
    return bonddims(bmp.M)
end

"""
    bonddims(bmp::BMP, i::Integer)

Returns the `i`-th bond dimension in the BMP.
"""
function bonddims(bmp::BMP, i::Integer)
    return bonddims(bmp.M, i)
end

function max_dim(bmp::BareBMP)
    return maximum(length(m.rows) for m in bmp.M[:,1])
end

function max_dim(bmp::BMP)
    return max_dim(bmp.M)
end

function volume(bmp::BareBMP, R::Vector{<:Integer}=[0,1])
    return sum(length(m.rows) for m in bmp[:,1]) + length(R)
end

"""
    volume(bmp::BMP)

Returns the sum of all bond dimensions, including the terminal vector length.
"""
function volume(bmp::BMP)
    return volume(bmp.M, bmp.R)
end

function evalfunc(bmp::BareBMP, x::BitArray, R::Vector{<:Integer}, order::Vector{<:Integer})::BitArray
    n = size(bmp, 1)
    m = length(bmp[1,1].rows)
    n_samps = div(length(x), size(x, 1))
    x_ = reshape(x, (n, n_samps))
    result = BitArray(undef, (m, n_samps))
    mat = fill(RSMInt(0), m)
    prod_mats = Vector{Vector{RSMInt}}(undef, n)
    for j=1:n_samps
        for i=1:n
            prod_mats[i] = bmp[i, x_[order[i], j]+1].rows
        end
        mat .= prod_mats[1]
        for i=2:n
            next_mat = prod_mats[i]
            for k=1:m
                mat[k] = next_mat[mat[k]]
            end
        end
        for k=1:m
            result[k,j] = R[mat[k]]
        end
    end
    shape = [i for i in size(x)]
    shape[1] = m
    return reshape(result, Tuple(shape))
end

"""
    evalfunc(bmp::BMP, x::BitArray)

Evaluates the BMP for the inputs given in `x`. The first dimension of `x` must
have dimension `length(bmp)`. The first dimension of the output will be the
number of output bits of `bmp`. The rest of the dimensions of `x` and the
output match.
"""
function evalfunc(bmp::BMP, x::BitArray)::BitArray
    return evalfunc(bmp.M, x, bmp.R, bmp.order)
end

function evalfunc(bmp::BMP, x::Array{<:Integer})
    f = evalfunc(bmp, x .!= 0)
    return Array{eltype(x)}(f)
end

function save_bmp(fpath::String, bmp::BMP)
    f = open(fpath, "w")
    n = length(bmp)
    println(f, "$n")
    for var_idx in bmp.order
        print(f, "$var_idx ")
    end
    println(f)
    dims = bonddims(bmp)
    for chi in dims
        print(f, "$chi ")
    end
    println(f)
    nR = length(bmp.R)
    println(f, "$nR")
    nL = bonddims(bmp, 1)
    for i=1:nL
        print(f, "1 ")
    end
    println(f, "")
    for val in bmp.R
        print(f, "$val ")
    end
    println(f, "")
    for i=1:n, j=1:2 # Note the order of the loops
        for val in bmp.M[i,j].rows
            print(f, "$val ")
        end
        println(f, "")
    end
    close(f)
end

function load_bmp(fpath)
    f = open(fpath)
    n = parse(UInt32, readline(f))
    order = parse.(UInt32, split(readline(f)))
    dims = parse.(UInt32, split(readline(f)))
    nR = parse(UInt32, readline(f))
    L = parse.(RSMInt, split(readline(f)))
    R = parse.(RSMInt, split(readline(f)))
    M = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        ncols = nR
        if i < n
            ncols = dims[i+1]
        end
        rows = parse.(RSMInt, split(readline(f)))
        M[i,j] = RowSwitchMatrix(rows, ncols)
    end
    return BMP(M, R, order)
end

