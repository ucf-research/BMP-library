include("RowSwitchMatrix.jl")

struct BMP
    M::Matrix{RowSwitchMatrix}
    R::Vector{RSMInt}
    order::Vector{UInt32}
    position::Vector{UInt32}
    function BMP(M::Matrix{RowSwitchMatrix}, R::Vector{<:Integer}, order::Vector{<:Integer})
        position = fill(0, length(order))
        # position is the inverse of the permutation given in order
        position[order] .= order
        return new(M, R, order, position)
    end
end

# Constructors and initializers
function BMP(val::Integer, n::Integer)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n-1
        mats[i,:] = [RowSwitchMatrix(1), RowSwitchMatrix(1)]
    end
    mats[n,1] = RowSwitchMatrix(RSMInt[val+1], 2)
    mats[n,2] = RowSwitchMatrix(RSMInt[val+1], 2)
    return BMP(mats, [0,1], collect(1:n))
end

function BMP(val::Integer, order::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n-1
        mats[i,:] = [RowSwitchMatrix(1), RowSwitchMatrix(1)]
    end
    mats[n,1] = RowSwitchMatrix(RSMInt[val+1], 2)
    mats[n,2] = RowSwitchMatrix(RSMInt[val+1], 2)
    return BMP(mats, [0,1], copy(order))
end

function BMP_bitline(xi::Integer, n::Integer)
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
    return BMP(mats, [0,1], collect(1:n))
end

function BMP_bitline(i::Integer, order::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    i = 1
    while order[i] != xi
        mats[i,:] = [RowSwitchMatrix(1), RowSwitchMatrix(1)]
        i += 1
    end
    mats[i,1] = RowSwitchMatrix(RSMInt[1], 2)
    mats[i,2] = RowSwitchMatrix(RSMInt[2], 2)
    for j=i+1:n
        mats[j,:] = [RowSwitchMatrix(2), RowSwitchMatrix(2)]
    end
    return BMP(mats, [0,1], copy(order))
end

# Convenience functions
function Base.length(bmp::BMP)
    return size(bmp.M, 1)
end

function BMP_dims(bmp::BMP)
    return [length(m) for m in bmp.M[:,1]]
end

function BMP_dims(bmp::BMP, i::Integer)
    return length(bmp.M[i,1])
end

function BMP_complexity(bmp::BMP)
    return maximum(length(m) for m in bmp.M[:,1])
end

function BMP_volume(bmp::BMP)
    return sum(length(m) for m in bmp.M[:,1]) + length(bmp.R)
end

# Evaluation
function eval(bmp::BMP, x::BitVector)::BitVector
    mat = bmp.M[1, x[bmp.order[1]]+1]
    for i=2:length(bmp)
        mat = RSM_mult(mat, bmp.M[i, x[bmp.order[i]]+1])
    end
    return bmp.R[mat.rows] .== 1
end

function eval(bmp::BMP, x::Vector{<:Integer})
    f = eval(bmp, x .!= 0)
    return Vector{Int64}(f)
end

# Cleaning
function BMP_clean1_lrstep(mats::Matrix{RowSwitchMatrix})
    # Unique elements of the vertically concatenated matrices on the left
    U = Dict{RSMInt, RSMInt}()
    L0 = [get!(U, i, length(U)+1) for i in mats[1,1].rows]
    L1 = [get!(U, i, length(U)+1) for i in mats[1,2].rows]
    result = Matrix{RowSwitchMatrix}(undef, (2,2))
    result[1,1] = RowSwitchMatrix(L0, length(U))
    result[1,2] = RowSwitchMatrix(L1, length(U))
    # Transformed matrices on the right
    R0 = fill(RSMInt(0), length(U))
    R1 = fill(RSMInt(0), length(U))
    for k in keys(U)
        ind = U[k]
        R0[ind] = mats[2,1].rows[k]
        R1[ind] = mats[2,2].rows[k]
    end
    result[2,1] = RowSwitchMatrix(R0, mats[2,1].ncols)
    result[2,2] = RowSwitchMatrix(R1, mats[2,2].ncols)
    return result
end

function BMP_clean1_rlstep(mats::Matrix{RowSwitchMatrix})
    # Store unique rows of the matrix pair on the right in a dictionary
    U = Dict{Tuple{RSMInt, RSMInt}, RSMInt}()
    for pair in zip(mats[2,1].rows, mats[2,2].rows)
        get!(U, pair, length(U)+1)
    end
    # Elements of the new matrices on the left
    result = Matrix{RowSwitchMatrix}(undef, (2,2))
    L0 = [get(U, (mats[2,1].rows[i], mats[2,2].rows[i]), 0) for i in mats[1,1].rows]
    result[1,1] = RowSwitchMatrix(L0, length(U))
    L1 = [get(U, (mats[2,1].rows[i], mats[2,2].rows[i]), 0) for i in mats[1,2].rows]
    result[1,2] = RowSwitchMatrix(L1, length(U))
    # Elements of the new matrices on the right
    R0 = fill(RSMInt(0), length(U))
    R1 = fill(RSMInt(0), length(U))
    for k in keys(U)
        ind = U[k]
        R0[ind] = k[1]
        R1[ind] = k[2]
    end
    result[2,1] = RowSwitchMatrix(R0, mats[2,1].ncols)
    result[2,2] = RowSwitchMatrix(R1, mats[2,2].ncols)
    return result
end

function BMP_clean1(bmp::BMP)
    n = length(bmp)
    M = copy(bmp.M)
    # Left-to-right sweep: eliminate unused rows
    for i=1:n-1
        new_pair = BMP_clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    # Right-to-left sweep: eliminate duplicate equivalent rows
    S_ = RowSwitchMatrix(bmp.R .+ 1, 2) # Trick to convert R vector to a row switching matrix
    M[n,:] .= [RSM_mult(M[n,1], S_), RSM_mult(M[n,2], S_)]
    for i=n-1:-1:1
        new_pair = BMP_clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return BMP(M, [0,1], copy(bmp.order))
end

function BMP_clean2_step(mats::Matrix{RowSwitchMatrix})
end

function BMP_clean2(bmp::BMP)
end

# Operations
function BMP_apply_noclean(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
    n = length(bmp1)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        mats[i,j] = RSM_kron(bmp1.M[i,j], bmp2.M[i,j])
    end
    return BMP(mats, copy(htab), copy(bmp1.order))
end

function BMP_apply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
    return BMP_clean1(BMP_apply_noclean(bmp1, bmp2, htab))
end

function BMP_apply_noclean(bmps::Vector{BMP}, htab::Vector{<:Integer})
    n = length(bmps[1])
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        m = bmps[1].M[i,j]
        for k=2:length(bmps)
            mats[i,j] = RSM_kron(m, bmps[k].M[i,j])
        end
    end
    return BMP(mats, copy(htab), copy(bmps[1].order))
end

function BMP_apply(bmps::Vector{BMP}, htab::Vector{<:Integer})
    return BMP_clean1(BMP_apply_noclean(bmps, htab))
end

function BMP_swap(bmp::BMP, i::Integer)
end
