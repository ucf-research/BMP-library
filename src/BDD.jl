const BDDNode = Tuple{UInt32, UInt32, UInt32}
struct BDD
    nodes::Vector{BDDNode}
    outputs::Vector{UInt32}
    order::Vector{UInt32}
    position::Vector{UInt32}
    function BDD(nodes, outputs, order)
        position = fill(UInt32(0), length(order))
        position[order] .= order
        return new(nodes, copy(outputs), copy(order), position)
    end
end

function Base.length(bdd::BDD)
    return length(bdd.order)
end

function BDD_reduce(bdd::BDD)
    n = length(bdd)
    V = length(bdd.nodes)
    sorted = sortperm(bdd.nodes, by=x->x[1], rev=true)
    ids = fill(UInt32(0), V)
    id_map = Dict{Tuple{UInt32, UInt32, UInt32}, UInt32}()
    for node_ind in sorted
        var_idx, lchild, rchild = bdd.nodes[node_ind]
        if var_idx == n+1
            lid, rid = 0, rchild
        else
            lid, rid = ids[lchild], ids[rchild]
        end
        if lid == rid && var_idx != n+1
            ids[node_ind] = lid
        else
            ids[node_ind] = get!(id_map, (var_idx, lid, rid), length(id_map)+1)
        end
    end
    newnodes = Vector{BDDNode}(undef, length(id_map))
    for (k,v) in pairs(id_map)
        newnodes[v] = k
    end
    newouts = [ids[i] for i in bdd.outputs]
    return BDD(newnodes, newouts, bdd.order)
end

function eval(bdd::BDD, x::BitArray)
    n = length(bdd)
    m = length(bdd.outputs)
    n_samps = div(length(x), size(x, 1))
    x_ = reshape(x, (n, n_samps))
    result = BitArray(undef, (m, n_samps))
    for j=1:n_samps
        for k=1:m
            node = bdd.nodes[bdd.outputs[k]]
            while node[1] != n+1
                next_ref = x_[node[1], j] == 0 ? node[2] : node[3]
                node = bdd.nodes[next_ref]
            end
            result[k, j] = node[3]
        end
    end
    shape = [i for i in size(x)]
    shape[1] = m
    return reshape(result, Tuple(shape))
end
