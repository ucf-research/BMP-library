include("RowSwitchMatrix.jl")

const BareBMP = Matrix{RowSwitchMatrix}

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
function BMP_bare(val::Integer, n::Integer)::BareBMP
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n-1
        mats[i,:] = [RowSwitchMatrix(1), RowSwitchMatrix(1)]
    end
    mats[n,1] = RowSwitchMatrix(RSMInt[val+1], 2)
    mats[n,2] = RowSwitchMatrix(RSMInt[val+1], 2)
    return mats
end

function BMP(val::Integer, n::Integer)
    mats = BMP_bare(val, n)
    return BMP(mats, [0,1], collect(1:n))
end

function BMP(val::Integer, order::Vector{<:Integer})
    mats = BMP_bare(val, length(order))
    return BMP(mats, [0,1], copy(order))
end

function BMP_bitline_bare(xi::Integer, n::Integer)::BareBMP
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

function BMP_bitline(xi::Integer, n::Integer)
    mats = BMP_bitline_bare(xi, n)
    return BMP(mats, [0,1], collect(1:n))
end

function BMP_bitline(xi::Integer, order::Vector{<:Integer})
    pos = findfirst(order .== xi)
    mats = BMP_bitline_bare(pos, length(order))
    return BMP(mats, [0,1], copy(order))
end

# Convenience functions
function Base.length(bmp::BMP)
    return size(bmp.M, 1)
end

function BMP_dims(bmp::BareBMP)
    return [length(m.rows) for m in bmp[:,1]]
end

function BMP_dims(bmp::BareBMP, i::Integer)
    return length(bmp[i,1].rows)
end

function BMP_dims(bmp::BMP)
    return BMP_dims(bmp.M)
end

function BMP_dims(bmp::BMP, i::Integer)
    return BMP_dims(bmp.M, i)
end

function BMP_complexity(bmp::BareBMP)
    return maximum(length(m.rows) for m in bmp.M[:,1])
end

function BMP_complexity(bmp::BMP)
    return BMP_complexity(bmp.M)
end

function BMP_volume(bmp::BareBMP, R::Vector{<:Integer}=[0,1])
    return sum(length(m.rows) for m in bmp[:,1]) + length(R)
end

function BMP_volume(bmp::BMP)
    return BMP_volume(bmp.M, bmp.R)
end

# Evaluation
function eval(bmp::BareBMP, x::BitArray, R::Vector{<:Integer}, order::Vector{<:Integer})::BitArray
    n = size(bmp, 1)
    m = length(bmp[1,1].rows)
    n_samps = div(length(x), size(x, 1))
    x_ = reshape(x, (n, n_samps))
    result = BitArray(undef, (m, n_samps))
    mat = fill(RSMInt(0), m)
    for j=1:n_samps
        mat .= bmp[1, x_[order[1],j] + 1].rows
        for i=2:n
            RSM_mult_inplace(mat, bmp[i, x_[order[i], j] + 1])
        end
        for k=1:m
            result[k,j] = R[mat[k]]
        end
    end
    shape = [i for i in size(x)]
    shape[1] = m
    return reshape(result, Tuple(shape))
end

function eval(bmp::BMP, x::BitArray)::BitArray
    return eval(bmp.M, x, bmp.R, bmp.order)
end

function eval(bmp::BMP, x::Array{<:Integer})
    f = eval(bmp, x .!= 0)
    return Array{eltype(x)}(f)
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

function BMP_clean1_rl(bmp::BareBMP, R::Vector{<:Integer})::BareBMP
    n = size(bmp, 1)
    M = copy(bmp)
    S_ = RowSwitchMatrix(R .+ 1, 2) # Trick to convert R vector to a row switching matrix
    M[n,1] = RSM_mult(M[n,1], S_)
    M[n,2] = RSM_mult(M[n,2], S_)
    for i=n-1:-1:1
        new_pair = BMP_clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

function BMP_clean1_rl(bmp::BMP)
    return BMP(BMP_clean1_rl(bmp.M, bmp.R), [0,1], copy(bmp.order))
end

function BMP_clean1_lr(bmp::BareBMP)::BareBMP
    n = size(bmp, 1)
    M = copy(bmp)
    for i=1:n-1
        new_pair = BMP_clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

function BMP_clean1_lr(bmp::BMP)
    return BMP(BMP_clean1_lr(bmp.M), copy(bmp.R), copy(bmp.order))
end

function BMP_clean1(bmp::BareBMP, R::Vector{<:Integer})::BareBMP
    n = size(bmp, 1)
    M = copy(bmp)
    # Left-to-right sweep: eliminate unused rows
    for i=1:n-1
        new_pair = BMP_clean1_lrstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    # Right-to-left sweep: eliminate duplicate equivalent rows
    S_ = RowSwitchMatrix(R .+ 1, 2) # Trick to convert R vector to a row switching matrix
    M[n,:] .= [RSM_mult(M[n,1], S_), RSM_mult(M[n,2], S_)]
    for i=n-1:-1:1
        new_pair = BMP_clean1_rlstep(M[i:i+1,:])
        M[i:i+1,:] .= new_pair
    end
    return M
end

function BMP_clean1(bmp::BMP)
    return BMP(BMP_clean1(bmp.M, bmp.R), [0,1], copy(bmp.order))
end

function BMP_clean2_step(mats::Matrix{RowSwitchMatrix})
end

function BMP_clean2(bmp::BMP)
end

# APPLY and helper routines
function BMP_apply_R(R1::Vector{<:Integer}, R2::Vector{<:Integer}, htab::Vector{<:Integer})
    R = [htab[2*i+j+1] for (j,i) in Iterators.product(R2, R1)]
    R = reshape(R, length(R))
    return R
end

function BMP_apply_R(Rs::Vector{<:Vector{<:Integer}}, htab::Vector{<:Integer})
    nrows = prod(length.(R for R in Rs))
    inds = fill(0, nrows)
    stride = nrows
    for R in Rs
        inds .*= 2
        mr = length(R)
        stride = div(stride, mr)
        cnt = div(nrows, stride)
        for i=0:cnt-1
            inds[i*stride+1:i*stride+stride] .+= R[i % mr + 1]
        end
    end
    return htab[inds .+ 1]
end

function BMP_apply_noclean(bmp1::BareBMP, bmp2::BareBMP)::BareBMP
    n = size(bmp1, 1)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        mats[i,j] = RSM_kron(bmp1[i,j], bmp2[i,j])
    end
    return mats
end

function BMP_apply(bmp1::BareBMP, bmp2::BareBMP, htab::Vector{<:Integer})
    # Assume canonical R
    return BMP_clean1(BMP_apply_noclean(bmp1, bmp2), htab)
end

function BMP_apply(bmp1::BareBMP, bmp2::BareBMP, R1::Vector{<:Integer}, R2::Vector{<:Integer}, htab::Vector{<:Integer})
    return BMP_clean1(BMP_apply_noclean(bmp1, bmp2), BMP_apply_R(R1, R2, htab))
end

function BMP_apply_noclean(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
    mats = BMP_apply_noclean(bmp1.M, bmp2.M)
    return BMP(mats, BMP_apply_R(bmp1.R, bmp2.R, htab), copy(bmp1.order))
end

function BMP_apply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})
    return BMP_clean1(BMP_apply_noclean(bmp1, bmp2, htab))
end

function BMP_apply_noclean(bmps::Vector{BareBMP})::BareBMP
    n = size(bmps[1], 1)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        mats[i,j] = RSM_kron([M[i,j] for M in bmps])
    end
    return mats
end

function BMP_apply(bmps::Vector{BareBMP}, htab::Vector{<:Integer})
    return BMP_clean1(BMP_apply_noclean(bmps), htab)
end

function BMP_apply(bmps::Vector{BareBMP}, Rs::Vector{<:Vector{<:Integer}}, htab::Vector{<:Integer})
    return BMP_clean1(BMP_apply_noclean(bmps), BMP_apply_R(Rs, htab))
end

function BMP_apply_noclean(bmps::Vector{BMP}, htab::Vector{<:Integer})
    mats = BMP_apply_noclean([bmp.M for bmp in bmps])
    return BMP(mats, BMP_apply_R([bmp.R for bmp in bmps], htab), copy(bmps[1].order))
end

function BMP_apply(bmps::Vector{BMP}, htab::Vector{<:Integer})
    return BMP_clean1(BMP_apply_noclean(bmps, htab))
end

# Alternative APPLY algorithm, faster in a lot of cases
function BMP_minapply(
    bmp1::BareBMP,
    bmp2::BareBMP
)
    n = size(bmp1,1)
    U = [(RSMInt(1),RSMInt(1))]
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i in axes(bmp1,1)
        A = Dict{Tuple{RSMInt, RSMInt}, RSMInt}()
        S0 = [get!(A, (bmp1[i,1].rows[i1], bmp2[i,1].rows[i2]), length(A)+1)
            for (i1, i2) in U]
        S1 = [get!(A, (bmp1[i,2].rows[i1], bmp2[i,2].rows[i2]), length(A)+1)
            for (i1, i2) in U]
        mats[i,1] = RowSwitchMatrix(S0, length(A))
        mats[i,2] = RowSwitchMatrix(S1, length(A))
        U = fill((RSMInt(0), RSMInt(0)), length(A))
        for k in keys(A)
            U[A[k]] = k
        end
    end
    return (mats, U)
end

function BMP_minapply_R(
    U::Vector{Tuple{RSMInt, RSMInt}},
    R1::Vector{<:Integer},
    R2::Vector{<:Integer},
    htab::Vector{<:Integer}
)
    R = fill(0, length(U))
    for (i,pair) in enumerate(U)
        i1, i2 = pair
        val = 2 * R1[i1] + R2[i2]
        R[i] = htab[val+1]
    end
    return R
end

function BMP_minapply(bmp1::BMP, bmp2::BMP, htab::Vector{<:Integer})::BMP
    M, U = BMP_minapply(bmp1.M, bmp2.M)
    return BMP(M, BMP_minapply_R(U, bmp1.R, bmp2.R, htab), copy(bmp1.order))
end

function BMP_minapply(bmps::Vector{BareBMP})
    n = size(bmps[1],1)
    k = length(bmps)
    U = fill(RSMInt(1), (k,1))
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:size(bmps[1],1)
        A = Dict{Vector{RSMInt}, RSMInt}()
        S0 = [get!(A, [bmps[j][i,1].rows[U[j,r]] for j=1:k], length(A)+1)
            for r in axes(U,2)]
        S1 = [get!(A, [bmps[j][i,2].rows[U[j,r]] for j=1:k], length(A)+1)
            for r in axes(U,2)]
        mats[i,1] = RowSwitchMatrix(S0, length(A))
        mats[i,2] = RowSwitchMatrix(S1, length(A))
        U = fill(RSMInt(1), (k,length(A)))
        for key in keys(A)
            U[:, A[key]] .= key
        end
    end
    return (mats, U)
end

function BMP_minapply_R(
    U::Matrix{RSMInt},
    Rs::Vector{<:Vector{<:Integer}},
    htab::Vector{<:Integer}
)
    R = fill(0, size(U,2))
    k = size(U,1)
    for i=1:size(U,2)
        val = 0
        for j=1:k
            val = 2*val + Rs[k][U[j,i]]
        end
        R[i] = htab[val+1]
    end
    return R
end

function BMP_minapply(bmps::Vector{BMP}, htab::Vector{<:Integer})::BMP
    M, U = BMP_minapply([bmp.M for bmp in bmps])
    R = BMP_minapply_R(U, [bmp.R for bmp in bmps], htab)
    return BMP_clean1_rl(BMP(M, R, copy(bmps[1].order)))
end

function BMP_multiapply_noclean(
    chip::BareBMP,
    bits::Vector{<:Vector{<:Integer}}
)
    n = size(chip, 1)
    U = [(RSMInt(k), Vector{RSMInt}(bs)) for (k,bs) in enumerate(bits)]
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        A = Dict{Tuple{RSMInt, Vector{RSMInt}}, RSMInt}()
        S0 = [get!(A, (k, [chip[i,1].rows[j] for j in bs]), length(A)+1)
            for (k,bs) in U]
        S1 = [get!(A, (k, [chip[i,2].rows[j] for j in bs]), length(A)+1)
            for (k,bs) in U]
        mats[i,1] = RowSwitchMatrix(S0, length(A))
        mats[i,2] = RowSwitchMatrix(S1, length(A))
        U = Vector{Tuple{RSMInt, Vector{RSMInt}}}(undef, length(A))
        for k in keys(A)
            U[A[k]] = k
        end
    end
    return (mats, U)
end

function BMP_multiapply(
    chip::BMP,
    bits::Vector{<:Vector{<:Integer}},
    tabs::Vector{<:Vector{<:Integer}}
)
    mats, U = BMP_multiapply_noclean(chip.M, bits)
    R = fill(RSMInt(0), length(U))
    for (i, row) in enumerate(U)
        k, bs = row
        tabval = 0
        for bval in bs
            tabval = 2 * tabval + bval - 1
        end
        R[i] = tabs[k][tabval + 1]
    end
    return BMP(BMP_clean1_rl(mats, R), [0,1], copy(chip.order))
end

function BMP_layerapply(
    bmp::BMP,
    bits::Vector{<:Vector{<:Integer}},
    tabs::Vector{<:Vector{<:Integer}}
)
    n = length(bmp)
    bits_ = Vector{Vector{UInt8}}(undef, n)
    tabs_ = Vector{Vector{UInt8}}(undef, n)
    for (gbits, gtab) in zip(bits, tabs)
        gsize = length(gbits)
        for (s, bline) in zip(gsize-1:-1:0, gbits)
            bits_[bline] = gbits
            tabs_[bline] = (gtab .>> s) .& 1
        end
    end
    return BMP_multiapply(bmp, bits_, tabs_)
end

# SWAP functions
function BMP_swap!(bmp::BareBMP, i::Integer)
    mats = bmp[i:i+1,:]
    chi = length(mats[1,1].rows)
    # Unique elements
    U = Dict{Tuple{RSMInt,RSMInt}, RSMInt}()
    sizehint!(U, 2 * chi)
    Ls = [Vector{RSMInt}(undef, chi), Vector{RSMInt}(undef, chi)]
    for j=1:2
        m2 = mats[2,j].rows
        for i in 1:chi
            i1 = mats[1,1].rows[i]
            i2 = mats[1,2].rows[i]
            k = (m2[i1], m2[i2])
            dv = RSMInt(length(U) + 1)
            Ls[j][i] = get!(U, k, dv)
        end
    end
    # Left matrices
    bmp[i,1] = RowSwitchMatrix(Ls[1], length(U))
    bmp[i,2] = RowSwitchMatrix(Ls[2], length(U))
    # Right matrices
    R0 = fill(RSMInt(0), length(U))
    R1 = fill(RSMInt(1), length(U))
    for (k,v) in pairs(U)
        R0[v] = k[1]
        R1[v] = k[2]
    end
    bmp[i+1,1] = RowSwitchMatrix(R0, mats[2,1].ncols)
    bmp[i+1,2] = RowSwitchMatrix(R1, mats[2,2].ncols)
end

function BMP_swap2!(bmp::BareBMP, i::Integer)
    mats = bmp[i:i+1,:]
    # Unique elements
    chi = length(mats[1,1].rows)
    U = Vector{Tuple{RSMInt, RSMInt}}(undef, 2 * chi)
    for j=1:2
        m2 = mats[2,j].rows
        for i=1:chi
            i1 = mats[1,1].rows[i]
            i2 = mats[1,2].rows[i]
            U[(j-1) * chi + i] = (m2[i1], m2[i2])
        end
    end
    unique!(sort!(U))
    # Left matrices
    L0 = [searchsortedfirst(U, (mats[2,1].rows[i1], mats[2,1].rows[i2]))
        for (i1,i2) in zip(mats[1,1].rows, mats[1,2].rows)]
    bmp[i,1] = RowSwitchMatrix(L0, length(U))
    L1 = [searchsortedfirst(U, (mats[2,2].rows[i1], mats[2,2].rows[i2]))
        for (i1,i2) in zip(mats[1,1].rows, mats[1,2].rows)]
    bmp[i,2] = RowSwitchMatrix(L1, length(U))
    # Right matrices
    R0 = fill(RSMInt(0), length(U))
    R1 = fill(RSMInt(1), length(U))
    for i in axes(U,1)
        R0[i] = U[i][1]
        R1[i] = U[i][2]
    end
    bmp[i+1,1] = RowSwitchMatrix(R0, mats[2,1].ncols)
    bmp[i+1,2] = RowSwitchMatrix(R1, mats[2,2].ncols)
end

function BMP_swap!(bmp::BMP, i::Integer)
    BMP_swap!(bmp.M, i)
    temp = bmp.order[i]
    bmp.order[i] = bmp.order[i+1]
    bmp.order[i+1] = temp
    bmp.position[bmp.order[i]] = i
    bmp.position[bmp.order[i+1]] = i+1
end

# JOIN
function BMP_join(bmps::Vector{BMP})
    n = size(bmps[1].M, 1)
    mats = Matrix{RowSwitchMatrix}(undef, (n, 2))
    for i=1:n, j=1:2
        mats[i,j] = RSM_join([bmp.M[i,j] for bmp in bmps])
    end
    R = fill(0, sum(length.(bmp.R for bmp in bmps)))
    stride = 0
    for bmp in bmps
        R[stride+1:stride+length(bmp.R)] .= bmp.R
        stride += length(bmp.R)
    end
    return BMP_clean1(BMP(mats, R, copy(bmps[1].order)))
end

# Save and load BMPs
function BMP_save(fpath::String, bmp::BMP)
    f = open(fpath, "w")
    n = length(bmp)
    println(f, "$n")
    for var_idx in bmp.order
        print(f, "$var_idx ")
    end
    println(f)
    dims = BMP_dims(bmp)
    for chi in dims
        print(f, "$chi ")
    end
    println(f)
    nR = length(bmp.R)
    println(f, "$nR")
    nL = BMP_dims(bmp, 1)
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

function BMP_load(fpath)
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
