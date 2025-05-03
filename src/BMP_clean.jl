function clean1_lrstep(mats::Matrix{RowSwitchMatrix})
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

function clean1_rlstep(mats::Matrix{RowSwitchMatrix})
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

function clean1_rl(bmp::BareBMP, R::Vector{<:Integer})::BareBMP
    n = size(bmp, 1)
    M = copy(bmp)
    S_ = RowSwitchMatrix(R .+ 1, 2) # Trick to convert R vector to a row switching matrix
    M[n,1] = mult(M[n,1], S_)
    M[n,2] = mult(M[n,2], S_)
    for i=n-1:-1:1
        new_pair = clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

"""
    clean1_rl(bmp::BMP)

Performs RTL-cleaning on `bmp`. The input BMP is not modified, the return value
is a new BMP.
"""
function clean1_rl(bmp::BMP)
    return BMP(clean1_rl(bmp.M, bmp.R), [0,1], copy(bmp.order))
end

function clean1_lr(bmp::BareBMP)::BareBMP
    n = size(bmp, 1)
    M = copy(bmp)
    for i=1:n-1
        new_pair = clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

"""
    clean1_lr(bmp::BMP)

Performs LTR-cleaning on `bmp`. The input BMP is not modified, the return value
is a new BMP.
"""
function clean1_lr(bmp::BMP)
    return BMP(clean1_lr(bmp.M), copy(bmp.R), copy(bmp.order))
end

function clean1(bmp::BareBMP, R::Vector{<:Integer})::BareBMP
    n = size(bmp, 1)
    M = copy(bmp)
    # Left-to-right sweep: eliminate unused rows
    for i=1:n-1
        new_pair = clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    # Right-to-left sweep: eliminate duplicate equivalent rows
    S_ = RowSwitchMatrix(R .+ 1, 2) # Trick to convert R vector to a row switching matrix
    M[n,:] .= [mult(M[n,1], S_), mult(M[n,2], S_)]
    for i=n-1:-1:1
        new_pair = clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

"""
    clean1(bmp::BMP)

Performs LTR-cleaning followed by RTL-cleaning on `bmp`. The input BMP is not
modified, the return value is a new BMP.
"""
function clean1(bmp::BMP)
    return BMP(clean1(bmp.M, bmp.R), [0,1], copy(bmp.order))
end
