const BDDSize = UInt32
const BDDVar = UInt16
struct BDDNode
    var_idx::BDDVar
    lchild::BDDSize
    hchild::BDDSize
end
struct BDD
    nodes::Vector{BDDNode}
    outputs::Vector{BDDSize}
    order::Vector{BDDVar}
    position::Vector{BDDVar}
    levels::Vector{Set{BDDSize}}
    ref_count::Vector{BDDSize}
    gc_refs::Set{BDDSize}
    function BDD(nodes, outputs, order)
        n = length(order)
        V = length(nodes)
        # Ensure the consistency of order and position
        position = fill(BDDVar(0), length(order))
        position[order] .= order
        # Ensure the consistency of levels and ref_count with nodes
        levels = [Set{BDDSize}() for _ in 1:length(order)+1]
        ref_count = [BDDSize(0) for _ in 1:length(nodes)]
        for (node_ref, node) in enumerate(nodes)
            if node.lchild > V || node.hchild > V
                error("BDD initialization failed: Non-existent nodes are referenced.")
            end
            if node.var_idx == n+1
                if node.lchild != 0 || (node.hchild != 0 && node.hchild != 1)
                    error("BDD initialization failed: Terminal nodes not properly filled.")
                end
            else
                lc_node = nodes[node.lchild]
                if lc_node.var_idx != n+1
                    if position[lc_node.var_idx] <= position[node.var_idx]
                        error("BDD initialization failed: Variable ordering not satisfied.")
                    end
                end
                hc_node = nodes[node.hchild]
                if hc_node.var_idx != n+1
                    if position[hc_node.var_idx] <= position[node.var_idx]
                        error("BDD initialization failed: Variable ordering not satisfied.")
                    end
                end
                ref_count[node.lchild] += 1
                ref_count[node.hchild] += 1
            end
            push!(levels[node.var_idx], node_ref)
        end
        for node_ref in outputs
            ref_count[node_ref] += 1
        end
        # Reuse empty indices
        gc_refs = Set{BDDSize}()
        for (node_ref, cnt) in enumerate(ref_count)
            if cnt == 0
                push!(gc_refs, node_ref)
            end
        end
        return new(nodes, outputs, order, position, levels, ref_count, gc_refs)
    end
end

function Base.length(bdd::BDD)
    return length(bdd.order)
end

function volume(bdd::BDD)
    return sum(cnt > 0 ? 1 : 0 for cnt in bdd.ref_count)
end

function gc_bdd!(bdd::BDD, var_idx::Integer)
    nodes = collect(bdd.levels[var_idx])
    for node_ref in nodes
        cnt = bdd.ref_count[node_ref]
        if cnt == 0
            delete!(bdd.levels[var_idx], node_ref)
            node = bdd.nodes[node_ref]
            bdd.ref_count[node.lchild] -= 1
            bdd.ref_count[node.hchild] -= 1
            push!(bdd.gc_refs, node_ref)
        end
    end
end

function gc_bdd!(bdd::BDD)
    for var_idx in axes(bdd.order,1)
        gc_bdd!(bdd, var_idx)
    end
end

function reduce_bdd!(bdd::BDD)
    n = length(bdd)
    ids = Dict(i => i for i in keys(bdd.nodes))
    id_map = Dict{BDDNode, BDDSize}()
    for var_idx in Iterators.flatten((n+1, Iterators.reverse(bdd.order)))
        for node_ref in bdd.levels[var_idx]
            node = bdd.nodes[node_ref]
            lid, hid = 0, node.hchild
            if var_idx != n+1
                lid, hid = ids[node.lchild], ids[node.hchild]
            end
            if var_idx != n+1 && lid == hid
                # passthrough node, remove
                ids[node_ref] = lid
            else
                newnode = BDDNode(node.var_idx, lid, hid)
                newid = get!(id_map, newnode, length(id_map)+1)
                ids[node_ref] = newid
            end
        end
    end
    V = length(id_map)
    resize!(bdd.nodes, V)
    resize!(bdd.ref_count, V)
    bdd.ref_count .= 0
    empty!(bdd.gc_refs)
    for level in bdd.levels
        empty!(level)
    end
    for (i,node_ref) in enumerate(bdd.outputs)
        bdd.outputs[i] = ids[node_ref]
        bdd.ref_count[bdd.outputs[i]] += 1
    end
    for (node, node_ref) in pairs(id_map)
        bdd.nodes[node_ref] = node
        push!(bdd.levels[node.var_idx], node_ref)
        if node.var_idx != n+1
            bdd.ref_count[node.lchild] += 1
            bdd.ref_count[node.hchild] += 1
        end
    end
    gc_bdd!(bdd)
end

function evalfunc(bdd::BDD, x::BitArray)
    n = length(bdd)
    m = length(bdd.outputs)
    n_samps = div(length(x), size(x, 1))
    x_ = reshape(x, (n, n_samps))
    result = BitArray(undef, (m, n_samps))
    for j=1:n_samps
        for k=1:m
            node = bdd.nodes[bdd.outputs[k]]
            while node.var_idx != n+1
                next_ref = x_[node.var_idx, j] == 0 ? node.lchild : node.hchild
                node = bdd.nodes[next_ref]
            end
            result[k, j] = node.hchild
        end
    end
    shape = [i for i in size(x)]
    shape[1] = m
    return reshape(result, Tuple(shape))
end

function create_node!(bdd::BDD)
    if length(bdd.gc_refs) > 0
        node_ref = first(bdd.gc_refs)
        delete!(bdd.gc_refs, node_ref)
        return node_ref
    end
    push!(bdd.nodes, BDDNode(0, 0, 0))
    push!(bdd.ref_count, 0)
    return length(bdd.nodes)
end

function create_node!(bdd::BDD, node::BDDNode)
    n = length(bdd)
    node_ref = create_node!(bdd)
    bdd.nodes[node_ref] = node
    if node.var_idx != n+1
        bdd.ref_count[node.lchild] += 1
        bdd.ref_count[node.hchild] += 1
    end
    push!(bdd.levels[node.var_idx], node_ref)
    return node_ref
end

function move_node!(bdd::BDD, node_ref::BDDSize, new_node::BDDNode)
    n = length(bdd)
    # Remove old version
    node = bdd.nodes[node_ref]
    if node.var_idx != n+1
        bdd.ref_count[node.lchild] -= 1
        bdd.ref_count[node.hchild] -= 1
    end
    delete!(bdd.levels[node.var_idx], node_ref)
    # Insert new version
    bdd.nodes[node_ref] = new_node
    if new_node.var_idx != n+1
        bdd.ref_count[new_node.lchild] += 1
        bdd.ref_count[new_node.hchild] += 1
    end
    push!(bdd.levels[new_node.var_idx], node_ref)
end

function swap!(bdd::BDD, i::Integer)
    x1 = bdd.order[i]
    x2 = bdd.order[i+1]
    nodes1 = collect(bdd.levels[x1])
    id_map = Dict(bdd.nodes[node_ref] => node_ref for node_ref in nodes1)
    for node_ref in nodes1
        node = bdd.nodes[node_ref]
        lc_ref, hc_ref = node.lchild, node.hchild
        lc, hc = bdd.nodes[lc_ref], bdd.nodes[hc_ref]
        if lc.var_idx == x2 && hc.var_idx == x2
            new_lc = BDDNode(x1, lc.lchild, hc.lchild)
            new_hc = BDDNode(x1, lc.hchild, hc.hchild)
        elseif lc.var_idx == x2 && hc.var_idx != x2
            new_lc = BDDNode(x1, lc.lchild, hc_ref)
            new_hc = BDDNode(x1, lc.hchild, hc_ref)
        elseif lc.var_idx != x2 && hc.var_idx == x2
            new_lc = BDDNode(x1, lc_ref, hc.lchild)
            new_hc = BDDNode(x1, lc_ref, hc.hchild)
        else
            continue
        end
        new_lc_ref = new_lc.lchild
        if new_lc.lchild != new_lc.hchild
            new_lc_ref = get(id_map, new_lc, 0)
            if new_lc_ref == 0
                new_lc_ref = create_node!(bdd, new_lc)
                id_map[new_lc] = new_lc_ref
            end
        end
        new_hc_ref = new_hc.lchild
        if new_hc.lchild != new_hc.hchild
            new_hc_ref = get(id_map, new_hc, 0)
            if new_hc_ref == 0
                new_hc_ref = create_node!(bdd, new_hc)
                id_map[new_hc] = new_hc_ref
            end
        end
        move_node!(bdd, node_ref, BDDNode(x2, new_lc_ref, new_hc_ref))
    end
    bdd.order[i] = x2
    bdd.position[x2] = i
    bdd.order[i+1] = x1
    bdd.position[x1] = i+1
    gc_bdd!(bdd, x2)
end

function save_bdd(fpath::String, bdd::BDD)
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
    for nd in bdd.nodes
        println(f, "$(nd.var_idx) $(nd.lchild) $(nd.hchild)")
    end
    close(f)
end

function load_bdd(fpath::String)
    f = open(fpath)
    n, n_out, n_nodes = parse.(BDDSize, split(readline(f)))
    order = parse.(BDDVar, split(readline(f)))
    outputs = parse.(BDDSize, split(readline(f)))
    nodes = Vector{BDDNode}(undef, n_nodes)
    for i=1:n_nodes
        var_idx, lchild, hchild = parse.(BDDSize, split(readline(f)))
        nodes[i] = BDDNode(var_idx, lchild, hchild)
    end
    return BDD(nodes, outputs, order)
end
