struct CustomSparseMatrix
    data::Vector{RSMInt}
    rbegin::Vector{UInt32}
    rend::Vector{UInt32}
end

struct CustomView
    data::Vector{RSMInt}
    rb::UInt32
    re::UInt32
end

function Base.hash(x::CustomView, h::UInt64)
    h = xor(0x77cfa1eef01bca90, h)
    # h = hash(x.re - x.rb, h)
    for i=x.re:-1:x.rb
        h = hash(x.data[i], h)
    end
    return h
end

function Base.isequal(x::CustomView, y::CustomView)
    if (x.re - x.rb) != (y.re - y.rb)
        return false
    end
    for (i1, i2) in zip(x.rb:x.re, y.rb:y.re)
        if x.data[i1] != y.data[i2]
            return false
        end
    end
    return true
end

function empty_sparse_mat(U::Dict{CustomView, RSMInt})
    chi = length(U)
    #
    rlen = fill(UInt32(0), chi)
    for (vw, k) in pairs(U)
        rlen[k] = vw.re - vw.rb + 1
    end
    shift = sum(rlen)
    #
    rbegin = fill(UInt32(0), 2 * chi)
    rend = fill(UInt32(0), 2 * chi)
    prev = 0
    for k=1:chi
        rbegin[k] = prev + 1
        rend[k] = prev + rlen[k]
        rbegin[k+chi] = rbegin[k] + shift
        rend[k+chi] = rend[k] + shift
        prev = rend[k]
    end
    newdata = fill(RSMInt(0), 2 * shift)
    return CustomSparseMatrix(newdata, rbegin, rend)
end

function propagate_mat(U::Dict{CustomView, RSMInt}, m::Matrix{RowSwitchMatrix})
    newmat = empty_sparse_mat(U)
    shift = div(length(newmat.data), 2)
    for (vw, k) in pairs(U)
        rdata = vw.data
        j = 0
        for i in newmat.rbegin[k]:newmat.rend[k]
            col = rdata[vw.rb + j]
            newmat.data[i] = m[j+1,1].rows[col]
            newmat.data[i+shift] = m[j+1,2].rows[col]
            j += 1
        end
    end
    return newmat
end

function propagate_mat(U::Dict{CustomView, RSMInt}, m::Vector{RowSwitchMatrix})
    newmat = empty_sparse_mat(U)
    shift = div(length(newmat.data), 2)
    for (vw, k) in pairs(U)
        rdata = vw.data
        j = 0
        for i in newmat.rbegin[k]:newmat.rend[k]
            newmat.data[i] = m[1].rows[rdata[vw.rb + j]]
            newmat.data[i+shift] = m[2].rows[rdata[vw.rb + j]]
            j += 1
        end
    end
    return newmat
end
