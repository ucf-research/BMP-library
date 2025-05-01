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

function swap!(bmp::BMP, i::Integer)
    swap!(bmp.M, i)
    temp = bmp.order[i]
    bmp.order[i] = bmp.order[i+1]
    bmp.order[i+1] = temp
    bmp.position[bmp.order[i]] = i
    bmp.position[bmp.order[i+1]] = i+1
end

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

