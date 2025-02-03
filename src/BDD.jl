const BDDNode = Tuple{UInt32, UInt32, UInt32}
struct BDD
    nodes::Vector{BDDNode}
    outputs::Vector{UInt32}
    order::Vector{UInt32}
    position::Vector{UInt32}
    function BDD(nodes, outputs, order)
        position = fill(UInt32(0), length(order))
        position[order] .= order
        return new(nodes, outputs, order, position)
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

function BDD_save(fpath::String, bdd::BDD)
    f = open(fpath, "w")
    n = length(bdd.order)
    n_out = length(bdd.outputs)
    n_nodes = length(bdd.nodes)
    println(f, "$n $n_out $n_nodes")
    for var_idx in bdd.order
        print(f, "$var_idx ")
    end
    println(f)
    for out_node in bdd.outputs
        print(f, "$out_node ")
    end
    println(f)
    for (var_idx, lchild, rchild) in bdd.nodes
        println(f, "$var_idx $lchild $rchild")
    end
    close(f)
end

function BDD_load(fpath::String)
    f = open(fpath)
    n, n_out, n_nodes = parse.(UInt32, split(readline(f)))
    order = parse.(UInt32, split(readline(f)))
    outputs = parse.(UInt32, split(readline(f)))
    nodes = Vector{BDDNode}(undef, n_nodes)
    for i=1:n_nodes
        var_idx, lchild, rchild = parse.(UInt32, split(readline(f)))
        nodes[i] = (var_idx, lchild, rchild)
    end
    return BDD(nodes, outputs, order)
end
