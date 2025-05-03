"""
    insert_var(bmp::BMP, p::Integer)

Returns a BMP that represents the same functions as `bmp`, but with an
additional (inert) input variable at position `p`. If `p` is not specified,
it is set to `length(bmp)+1`, i.e. the end of the matrix train. The label of
the new variable is always `length(bmp)+1`.
"""
function insert_var(bmp::BMP, p::Integer)
    n = length(bmp)
    M = Matrix{RowSwitchMatrix}(undef, (n+1,2))
    d = length(bmp.R)
    if p < n+1
        d = length(bmp.M[p,1].rows)
    end
    M[1:p-1,:] .= bmp.M[1:p-1,:]
    M[p+1:n+1,:] .= bmp.M[p:n,:]
    M[p,1] = RowSwitchMatrix(d)
    M[p,2] = RowSwitchMatrix(d)
    order = copy(bmp.order)
    insert!(order, p, n+1)
    return BMP(M, copy(bmp.R), order)
end

function insert_var(bmp::BMP)
    n = length(bmp)
    return insert_var(bmp, n+1)
end

"""
    restrict(bmp::BMP, var::Integer, val::Integer)

Return the BMP of the function that is obtained upon fixing the value of
variable `var` to `val` in the function represented by `bmp`. The variable with
the label `var` is not removed.
"""
function restrict(bmp::BMP, var::Integer, val::Integer)
    M = copy(bmp.M)
    p = bmp.position[var]
    mat = M[p,val+1]
    M[p,2-val] = RowSwitchMatrix(copy(mat.rows), mat.ncols)
    return BMP(clean1(M, bmp.R), [0,1], copy(bmp.order))
end

"""
    erase_var(bmp::BMP, var::Integer, val::Integer)

Similar to [`restrict`](@ref), but removes the inert variable from the BMP.
This causes a relabeling of all variables with labels greater than `var`, so
the labels in the resulting BMP are `1` through `length(bmp)-1`.
"""
function erase_var(bmp::BMP, var::Integer, val::Integer)
    n = length(bmp)
    p = bmp.position[var]
    M = vcat(bmp.M[1:p-1,:], bmp.M[p+1:n,:])
    mat = bmp.M[p,val+1]
    if p == n
        M[n-1,1] = mult(bmp.M[n-1,1], mat)
        M[n-1,2] = mult(bmp.M[n-1,2], mat)
    else
        M[p,1] = mult(mat, bmp.M[p+1,1])
        M[p,2] = mult(mat, bmp.M[p+1,2])
    end
    order = vcat(bmp.order[1:p-1], bmp.order[p+1:n])
    for (i,vl) in enumerate(order)
        if vl > var
            order[i] -= 1
        end
    end
    return BMP(clean1(M, bmp.R), [0,1], order)
end

function compose(bmp::BMP, var::Integer, sub::BMP)
    t1 = apply(sub, restrict(bmp, var, 0), [0,1,0,0])
    t2 = apply(sub, restrict(bmp, var, 1), [0,0,0,1])
    return apply(t1, t2, [0,1,1,1])
end

function swap!(bmp::BareBMP, i::Integer)
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

"""
    swap!(bmp::BMP, i::Integer)

Swaps the input variables at positions `i` and `i+1`. Note that unlike most
BMP functions, `swap!` modifies `bmp`.
"""
function swap!(bmp::BMP, i::Integer)
    swap!(bmp.M, i)
    temp = bmp.order[i]
    bmp.order[i] = bmp.order[i+1]
    bmp.order[i+1] = temp
    bmp.position[bmp.order[i]] = i
    bmp.position[bmp.order[i+1]] = i+1
end

"""
    reorder!(bmp::BMP, pord::Vector{<:Integer})

Changes the variable ordering of `bmp` to `pord` using local swaps via
[`swap!`](@ref).
"""
function reorder!(bmp, pord::Vector{<:Integer})
    for (dest, var) in enumerate(pord)
        src = bmp.position[var]
        for i=src-1:-1:dest
            swap!(bmp, i)
        end
    end
end

"""
    joinfuncs(bmps::Vector{BMP})

Generates a joint BMP of the functions represented by the elements of `bmps`.
The outputs bits of the BMPs in `bmps` are stacked on top of each other in the
order they are given in `bmps`.
"""
function joinfuncs(bmps::Vector{BMP})
    n = size(bmps[1].M, 1)
    mats = Matrix{RowSwitchMatrix}(undef, (n, 2))
    for i=1:n, j=1:2
        mats[i,j] = dsum([bmp.M[i,j] for bmp in bmps])
    end
    R = fill(0, sum(length.(bmp.R for bmp in bmps)))
    stride = 0
    for bmp in bmps
        R[stride+1:stride+length(bmp.R)] .= bmp.R
        stride += length(bmp.R)
    end
    return clean1(BMP(mats, R, copy(bmps[1].order)))
end
