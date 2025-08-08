# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

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
    for (k, ind) in pairs(U)
        R0[ind] = mats[2,1].rows[k]
        R1[ind] = mats[2,2].rows[k]
    end
    result[2,1] = RowSwitchMatrix(R0, mats[2,1].ncols)
    result[2,2] = RowSwitchMatrix(R1, mats[2,2].ncols)
    return result
end

function clean1_rlstep(mats::Matrix{RowSwitchMatrix})
    # SU decomposition: S as an array and U as a dictionary
    U = Dict{Tuple{RSMInt, RSMInt}, RSMInt}()
    S = Vector{RSMInt}(undef, length(mats[2,1].rows))
    for (i, pair) in enumerate(zip(mats[2,1].rows, mats[2,2].rows))
        S[i] = get!(U, pair, length(U)+1)
    end
    # Elements of the new matrices on the left
    result = Matrix{RowSwitchMatrix}(undef, (2,2))
    L0 = [S[i] for i in mats[1,1].rows]
    result[1,1] = RowSwitchMatrix(L0, length(U))
    L1 = [S[i] for i in mats[1,2].rows]
    result[1,2] = RowSwitchMatrix(L1, length(U))
    # NOTE: L0, L1 can be obtained directly from U with no S, but using S
    # turns out to be more efficient
    #
    # Elements of the new matrices on the right
    R0 = fill(RSMInt(0), length(U))
    R1 = fill(RSMInt(0), length(U))
    for (k, ind) in pairs(U)
        R0[ind] = k[1]
        R1[ind] = k[2]
    end
    result[2,1] = RowSwitchMatrix(R0, mats[2,1].ncols)
    result[2,2] = RowSwitchMatrix(R1, mats[2,2].ncols)
    return result
end

function clean1_rl(bmp::BareBMP, R::AbstractArray)
    n = size(bmp, 1)
    M = copy(bmp)
    S0_rows = Vector{RSMInt}(R)
    S0_rows .+= 1 # Trick to convert R vector to a row switching matrix
    S0 = RowSwitchMatrix(S0_rows, 2)
    M[n,1] = mult(M[n,1], S0)
    M[n,2] = mult(M[n,2], S0)
    for i=n-1:-1:1
        new_pair = clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

"""
    clean1_rl(bmp::BMP)

Performs RTL-cleaning on `bmp`. `bmp` is not modified, the return value is a
new BMP.
"""
function clean1_rl(bmp::BMP)
    return BMP(clean1_rl(bmp.M, bmp.R), RSMInt[0,1], copy(bmp.order))
end

function clean1_lr(bmp::BareBMP)
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

Performs LTR-cleaning on `bmp`. `bmp` is not modified, the return value is a
new BMP.
"""
function clean1_lr(bmp::BMP)
    return BMP(clean1_lr(bmp.M), copy(bmp.R), copy(bmp.order))
end

function clean1(bmp::BareBMP, R::Vector{<:Integer})
    n = size(bmp, 1)
    M = copy(bmp)
    # Left-to-right sweep: eliminate unused rows
    for i=1:n-1
        new_pair = clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    # Right-to-left sweep: eliminate duplicate equivalent rows
    S0_rows = Vector{RSMInt}(R)
    S0_rows .+= 1 # Trick to convert R vector to a row switching matrix
    S0 = RowSwitchMatrix(S0_rows, 2)
    M[n,1] = mult(M[n,1], S0)
    M[n,2] = mult(M[n,2], S0)
    for i=n-1:-1:1
        new_pair = clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

"""
    clean1(bmp::BMP)

Performs LTR-cleaning followed by RTL-cleaning on `bmp`. `bmp` is not
modified, instead, a new BMP is returned.
"""
function clean1(bmp::BMP)
    return BMP(clean1(bmp.M, bmp.R), RSMInt[0,1], copy(bmp.order))
end
