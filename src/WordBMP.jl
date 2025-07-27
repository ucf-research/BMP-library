# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

struct WordBMP
    M::Matrix{RowSwitchMatrix}
    R::Vector{RSMInt}
    order::Matrix{Int16}
    wpos::Vector{Int16}
    bpos::Vector{Int16}
    isize::Int16
    osize::Int16
    function WordBMP(
        M::Matrix{RowSwitchMatrix},
        R::Vector{<:Integer},
        order::Matrix{<:Integer},
        isize::Integer,
        osize::Integer
    )
        # Enforce correct dimensions & wpos, bpos entries
        if 2^isize != size(M,2)
            throw(AssertionError("Input word size does not match the number of matrices per site."))
        end
        n = length(order)
        wpos = fill(Int16(0), n)
        bpos = fill(Int16(0), n)
        for i in axes(order,1), j in axes(order,2)
            b = order[i, j]
            bpos[b] = i
            wpos[b] = j
        end
        return new(M, R, order, wpos, bpos, isize, osize)
    end
end

function consolidate_outputs(
    M::BareBMP,
    k::Integer
)
    n = size(M, 1)
    m = length(M[1,1].rows)
    if m % k != 0
        throw(AssertionError("Word size is not a divisor of the number of outputs."))
    end
    U = [collect(i:i+k-1) for i=1:k:m]
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        A = Dict{Vector{RSMInt}, RSMInt}()
        S0 = [get!(A, [M[i,1].rows[j] for j in bs], length(A)+1)
            for bs in U]
        S1 = [get!(A, [M[i,2].rows[j] for j in bs], length(A)+1)
            for bs in U]
        mats[i,1] = RowSwitchMatrix(S0, length(A))
        mats[i,2] = RowSwitchMatrix(S1, length(A))
        U = Vector{Vector{RSMInt}}(undef, length(A))
        for (k,v) in pairs(A)
            U[v] = k
        end
    end
    return (mats, U)
end

function consolidate_term(U::Vector{<:Vector{<:Integer}}, R::Vector{<:Integer})
    new_R = Vector{RSMInt}(undef, length(U))
    for (i,bs) in enumerate(U)
        val = 0
        for bind in bs
            val = 2 * val + R[bind]
        end
        new_R[i] = val
    end
    return new_R
end

function consolidate_inputs(
    M::Matrix{RowSwitchMatrix},
    k::Integer
)
    n = size(M, 1)
    if n % k != 0
        throw(AssertionError("Word size is not a divisor of the number of inputs."))
    end
    n_words = div(n, k)
    result = Matrix{RowSwitchMatrix}(undef, (n_words, 2^k))
    for wi=1:n_words
        wstart = (wi - 1) * k
        for val=0:2^k-1
            bval = val >> (k-1) & 1
            rsm = M[wstart + 1, bval + 1]
            mat = copy(rsm.rows)
            ncols = rsm.ncols
            for bi=2:k
                bval = val >> (k - bi) & 1
                rsm = M[wstart + bi, bval + 1]
                mult_inplace(mat, rsm)
                ncols = rsm.ncols
            end
            result[wi, val+1] = RowSwitchMatrix(mat, ncols)
        end
    end
    return result
end

function WordBMP(bmp::BMP, isize::Integer, osize::Integer)
    mats_, U = consolidate_outputs(bmp.M, osize)
    mats = consolidate_inputs(mats_, isize)
    R = consolidate_term(U, bmp.R)
    n = length(bmp.order)
    order = reshape(copy(bmp.order), (isize, div(n, isize)))
    return WordBMP(mats, R, order, isize, osize)
end

function break_inputs(M::Matrix{RowSwitchMatrix}, k::Integer)
    n_words = size(M, 1)
    bmp_mats = Matrix{RowSwitchMatrix}(undef, (n_words * k, 2))
    for wi=1:n_words
        n_rows = length(M[wi,1].rows)
        for bi=1:k-1
            bmp_mats[(wi-1)*k+bi,1] = RowSwitchMatrix(collect(1:n_rows), 2*n_rows)
            bmp_mats[(wi-1)*k+bi,2] = RowSwitchMatrix(collect(n_rows+1:2*n_rows), 2*n_rows)
            n_rows *= 2
        end
        blocks0 = Vector{Vector{RSMInt}}(undef, 2^(k-1))
        blocks1 = Vector{Vector{RSMInt}}(undef, 2^(k-1))
        for val=0:2^(k-1)-1
            ind = 0
            for i=0:k-1 # val has k-1 bits, we reverse k bits
                ind |= (val >> i & 1) << (k-1-i)
            end
            blocks0[val+1] = M[wi, ind+1].rows
            blocks1[val+1] = M[wi, ind+2].rows
        end
        ncols = M[wi, 1].ncols
        bmp_mats[(wi-1)*k+k, 1] = RowSwitchMatrix(collect(Iterators.flatten(blocks0)), ncols)
        bmp_mats[(wi-1)*k+k, 2] = RowSwitchMatrix(collect(Iterators.flatten(blocks1)), ncols)
    end
    return bmp_mats
end

function break_outputs(R::Vector{<:Integer}, k::Integer, n_words::Integer)
    new_R = reshape([w >> i & 1 for w in R, i=k-1:-1:0], (k * length(R)))
    n_outputs = k * n_words
    U = collect(Iterators.flatten([i:n_words:n_outputs for i=1:n_words]))
    return (U, new_R)
end

function BMP(wbmp::WordBMP)
    n_words = length(wbmp.M[1,1].rows)
    M_ = break_inputs(wbmp.M, wbmp.isize)
    M = [dsum([m for k=1:wbmp.osize]) for m in M_]
    U, R = break_outputs(wbmp.R, wbmp.osize, n_words)
    M[1,1] = RowSwitchMatrix(M[1,1].rows[U], M[1,1].ncols)
    M[1,2] = RowSwitchMatrix(M[1,2].rows[U], M[1,2].ncols)
    return clean1(BMP(M, R, reshape(wbmp.order, size(M,1))))
end

function bonddims(bmp::WordBMP) 
    return length.(m.rows for m in bmp.M[:,1])
end

function matrix_volume(bmp::WordBMP)
    return sum(bonddims(bmp))
end

function scaled_volume(bmp::WordBMP)
    return 2^bmp.osize * matrix_volume(bmp)
end

function volume(bmp::WordBMP)
    return matrix_volume(bmp) + length(bmp.R)
end

function compute_words(k::Integer, nw::Integer, xs::BitVector, order::Matrix{<:Integer})
    result = fill(0, nw)
    for i=1:nw
        for bi=1:k
            result[i] *= 2
            result[i] += xs[order[bi, i]]
        end
    end
    return result
end

function evalfunc(bmp::WordBMP, x::BitArray)
    n_words = size(bmp.M, 1)
    kin = bmp.isize
    kout = bmp.osize
    m_words = length(bmp.M[1,1].rows)
    n_samps = div(length(x), size(x, 1))
    #
    x_ = reshape(x, (n_words * kin, n_samps))
    result = BitArray(undef, (m_words * kout, n_samps))
    mat = fill(RSMInt(0), m_words)
    for j=1:n_samps
        sample_words = compute_words(kin, n_words, x_[:,j], bmp.order)
        mat .= bmp.M[1, sample_words[1] + 1].rows
        for i=2:n_words
            mult_inplace(mat, bmp.M[i, sample_words[i] + 1])
        end
        for i=1:m_words, bi=1:kout
            result[(i-1)*kout + bi, j] = bmp.R[mat[i]] >> (kout - bi) & 1
        end
    end
    shape = [i for i in size(x)]
    shape[1] = kout * m_words
    return reshape(result, Tuple(shape))
end

function word_clean1_rlstep(mats::Matrix{RowSwitchMatrix})
    K = size(mats, 2)
    # Store unique rows of the matrix pair on the right in a dictionary
    chi = length(mats[2,1].rows)
    U = Dict{Vector{RSMInt}, RSMInt}()
    sizehint!(U, chi)
    rmap = fill(0, chi)
    for r=1:chi
        rmap[r] = get!(U, [mats[2,i].rows[r] for i=1:K], length(U)+1)
    end
    # Elements of the new matrices on the left
    result = Matrix{RowSwitchMatrix}(undef, (2,K))
    for li=1:K
        L = fill(RSMInt(0), length(mats[1,li].rows))
        for (ridx,i) in enumerate(mats[1,li].rows)
            L[ridx] = rmap[i]
        end
        # L = [get(U, [mats[2,ri].rows[i] for ri=1:K], 0) for i in mats[1,li].rows]
        result[1,li] = RowSwitchMatrix(L, length(U))
    end
    # Elements of the new matrices on the right
    Rs = [fill(RSMInt(0), length(U)) for _ in 1:K]
    for (kd,v) in pairs(U)
        for i=1:K
            Rs[i][v] = kd[i]
        end
    end
    for i=1:K
        result[2,i] = RowSwitchMatrix(Rs[i], mats[2,i].ncols)
    end
    return result
end

function word_clean1_rl(mats::Matrix{RowSwitchMatrix}, R::Vector{<:Integer})
    n = size(mats, 1)
    K = size(mats, 2)
    M = copy(mats)
    S_ = RowSwitchMatrix(R .+ 1, K) # Trick to convert R vector to a row switching matrix
    for j=1:K
        M[n,j] = mult(M[n,j], S_)
    end
    for i=n-1:-1:1
        new_pair = word_clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

function clean1_rl(bmp::WordBMP)
    n = size(bmp.M, 1)
    K = 2^bmp.osize
    M = copy(bmp.M)
    S_ = RowSwitchMatrix(bmp.R .+ 1, K) # Trick to convert R vector to a row switching matrix
    for j=1:K
        M[n,j] = mult(M[n,j], S_)
    end
    for i=n-1:-1:1
        new_pair = word_clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return WordBMP(M, collect(0:K-1), copy(bmp.order), bmp.isize, bmp.osize)
end

function word_clean1_lrstep(mats::Matrix{RowSwitchMatrix})
    K = size(mats, 2)
    # Unique elements of the vertically concatenated matrices on the left
    chi = length(mats[1,1].rows)
    U = Dict{RSMInt, RSMInt}()
    sizehint!(U, 2*chi)
    Ls = Vector{Vector{RSMInt}}(undef, K)
    for j=1:K
        Ls[j] = [get!(U, i, length(U)+1) for i in mats[1,j].rows]
    end
    result = Matrix{RowSwitchMatrix}(undef, (2,K))
    for j=1:K
        result[1,j] = RowSwitchMatrix(Ls[j], length(U))
    end
    # Transformed matrices on the right
    Rs = [fill(RSMInt(0), length(U)) for _ in 1:K]
    for (k,ind) in pairs(U)
        for j=1:K
            Rs[j][ind] = mats[2,j].rows[k]
        end
    end
    for j=1:K
        result[2,j] = RowSwitchMatrix(Rs[j], mats[2,j].ncols)
    end
    return result
end

function word_clean1_lr(mats::Matrix{RowSwitchMatrix})
    n = size(mats, 1)
    M = copy(mats)
    for i=1:n-1
        new_pair = word_clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

function clean1_lr(bmp::WordBMP)
    n = size(bmp, 1)
    M = copy(bmp.M)
    for i=1:n-1
        new_pair = word_clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return WordBMP(M, copy(bmp.R), copy(bmp.order), bmp.isize, bmp.osize)
end

function word_clean1(mats::Matrix{RowSwitchMatrix}, R::Vector{<:Integer})
    n = size(mats, 1)
    K = size(mats, 2)
    M = copy(mats)
    # Left-to-right sweep: eliminate unused rows
    for i=1:n-1
        new_pair = word_clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    # Right-to-left sweep: eliminate duplicate equivalent rows
    S_ = RowSwitchMatrix(R .+ 1, K) # Trick to convert R vector to a row switching matrix
    for j=1:K
        M[n,j] = mult(M[n,j], S_)
    end
    for i=n-1:-1:1
        new_pair = word_clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

function clean1(bmp::WordBMP)
    n = size(bmp.M, 1)
    K = size(bmp.M, 2)
    M = copy(bmp.M)
    # Left-to-right sweep: eliminate unused rows
    for i=1:n-1
        new_pair = word_clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    # Right-to-left sweep: eliminate duplicate equivalent rows
    S_ = RowSwitchMatrix(bmp.R .+ 1, K) # Trick to convert R vector to a row switching matrix
    for j=1:K
        M[n,j] = mult(M[n,j], S_)
    end
    for i=n-1:-1:1
        new_pair = word_clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return WordBMP(M, collect(0:K-1), copy(bmp.order), bmp.isize, bmp.osize)
end

function joinfuncs(bmps::Vector{<:WordBMP})
    n = size(bmps[1].M, 1)
    K = size(bmps[1].M, 2)
    M = Matrix{RowSwitchMatrix}(undef, (n, K))
    join_mats = Vector{RowSwitchMatrix}(undef, length(bmps))
    for i=1:n, j=1:K
        for (mi, bmp) in enumerate(bmps)
            join_mats[mi] = bmp.M[i, j]
        end
        M[i,j] = dsum(join_mats)
    end
    nR = sum(length.(bmp.R for bmp in bmps))
    R = fill(RSMInt(0), nR)
    offset = 0
    for bmp in bmps
        sz = length(bmp.R)
        R[offset+1:offset+sz] .= bmp.R
        offset += sz
    end
    M1 = word_clean1_lr(M)
    M2 = word_clean1_rl(M1, R)
    return WordBMP(M2, collect(0:K-1), copy(bmps[1].order), bmps[1].isize, bmps[1].osize)
end

function canonize_terminal!(bmp::WordBMP)
    n = size(bmp.M, 1)
    K = size(bmp.M, 2)
    S = RowSwitchMatrix(bmp.R .+ 1, K) # Trick to convert R vector to a row switching matrix
    for j=1:K
        bmp.M[n,j] = mult(bmp.M[n,j], S)
    end
    resize!(bmp.R, K)
    for j=1:K
        bmp.R[j] = j-1
    end
end

function transform_outputs!(bmp::WordBMP, perm::Vector{<:Integer})
    for (i, val) in enumerate(bmp.R)
        bmp.R[i] = perm[val+1]
    end
    canonize_terminal!(bmp)
end

function transform_inputs!(bmp::WordBMP, a::Integer, perm::Vector{<:Integer})
    M = bmp.M[a,:]
    K = size(bmp.M, 2)
    for val=0:K-1
        bmp.M[a, val+1] = M[perm[val+1]+1]
    end
end
