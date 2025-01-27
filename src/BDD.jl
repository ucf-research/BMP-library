const BDDNode = Tuple{UInt32, UInt32, UInt32}
struct BDD
    nodes::Vector{BDDNode}
    order::Vector{UInt32}
    position::Vector{UInt32}
    function BDD(nodes, order)
        position = fill(UInt32(0), length(order))
        position[order] .= order
        return new(nodes, order, position)
    end
end

function Base.length(bdd::BDD)
    return length(bdd.order)
end

function BDD_reduce!(bdd::BDD)
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
    resize!(bdd.nodes, length(id_map))
    for (k,v) in pairs(id_map)
        bdd.nodes[v] = k
    end
end
