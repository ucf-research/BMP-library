function get_swap_position(n::Integer, i::Integer)
    if n == 2
        return 1
    end
    l = i % n
    k = div(i, n)
    if l != 0 && k % 2 == 0
        return n-l
    elseif l != 0 && k % 2 == 1
        return l
    elseif k % 2 == 1
        return get_swap_position(n-1, k) + 1
    else
        return get_swap_position(n-1, k)
    end
end

"""
    brute_force!(bmp::BMP)

Optimizes the variable ordering of `bmp` by brute-force iterating over all
possible permutations of variables via local swaps. This is a very inefficient
method that mainly exists to check other optimization methods.
"""
function brute_force!(bmp::BMP)
    min_order = copy(bmp.order)
    min_vol = volume(bmp)
    n = length(bmp)
    N = prod(1:n)
    for i=1:N-1
        p = get_swap_position(n, i)
        swap!(bmp, p)
        vol = volume(bmp)
        if vol < min_vol
            min_order .= bmp.order
            min_vol = vol
        end
    end
    reorder!(bmp, min_order)
end

# A heap based priority queue used for exact minimization
struct CustomHeap
    costs::Vector{Tuple{UInt64, BitVector}}
    states::Dict{BitVector, UInt64}
end

function CustomHeap()
    return CustomHeap(Tuple{UInt64, BitVector}[], Dict{BitVector, UInt64}())
end

function switch_elems!(heap::CustomHeap, i::Integer, j::Integer)
    c1, s1 = heap.costs[i]
    c2, s2 = heap.costs[j]
    heap.states[s1] = j
    heap.states[s2] = i
    heap.costs[i] = (c2, s2)
    heap.costs[j] = (c1, s1)
end

function move_down!(heap::CustomHeap, i::Integer)
    n_elems = length(heap.costs)
    cost_i = heap.costs[i][1]
    ci1, ci2 = 2*i, 2*i+1
    while ci1 <= n_elems
        if ci2 <= n_elems
            ci = heap.costs[ci1][1] < heap.costs[ci2][1] ? ci1 : ci2
        else
            ci = ci1
        end
        if heap.costs[ci][1] < cost_i
            switch_elems!(heap, i, ci)
            i = ci
        else
            break
        end
        ci1, ci2 = 2*i, 2*i + 1
    end
end

function move_up!(heap::CustomHeap, i::Integer)
    pi = div(i, 2)
    cost_i = heap.costs[i][1]
    while pi > 0
        if heap.costs[pi][1] >= cost_i
            switch_elems!(heap, pi, i)
            i = pi
            pi = div(i, 2)
        else
            break
        end
    end
end

function Base.haskey(heap::CustomHeap, state::BitVector)
    return haskey(heap.states, state)
end

function update!(heap::CustomHeap, state::BitVector, cost::Integer)
    if haskey(heap.states, state)
        ind = heap.states[state]
        pcost, state_ = heap.costs[ind]
        heap.costs[ind] = (cost, state_)
        if pcost > cost
            move_up!(heap, ind)
        else
            move_down!(heap, ind)
        end
    else
        state_ = copy(state)
        push!(heap.costs, (cost, state_))
        heap.states[state_] = length(heap.costs)
        move_up!(heap, length(heap.costs))
    end
end

function Base.pop!(heap::CustomHeap)
    if length(heap.costs) == 0
        throw(ArgumentError("Heap must be non-empty."))
    end
    c, s = heap.costs[1]
    delete!(heap.states, s)
    head = pop!(heap.costs)
    if length(heap.costs) > 0
        heap.costs[1] = head
        heap.states[head[2]] = 1
        move_down!(heap, 1)
    end
    return (c,s)
end

# Other helper functions for exact minimization
function partial_order!(bmp::BMP, pord::Vector{<:Integer})
    for (dest, var) in enumerate(pord)
        src = bmp.position[var]
        for i=src-1:-1:dest
            swap!(bmp, i)
        end
    end
end

function reconstruct_ordering!(
    bmp::BMP,
    pord::Vector{<:Integer},
    min_vol::Integer,
    min_order::Vector{<:Integer}
)
    for (dest, var) in enumerate(pord)
        src = bmp.position[var]
        for i=src-1:-1:dest
            swap!(bmp, i)
            vol = volume(bmp)
            if vol < min_vol
                min_vol = vol
                min_order .= bmp.order
            end
        end
    end
    return min_vol
end

function trace_path(state::BitVector, prev::Dict{BitVector, UInt64})
    state_ = copy(state)
    k = count(state_)
    pord = fill(0, k)
    for i=k:-1:1
        var = prev[state_]
        pord[i] = var
        state_[var] = false
    end
    return pord
end

function compute_heuristic(bdim::Integer, nr::Integer)
    total = 0
    prev_bond = bdim
    for _ in 1:nr
        d = 1
        while d * d < prev_bond
            d += 1
        end
        total += d
        prev_bond = d
    end
    return total
end

function basic_exact_minimize!(bmp::BMP)
    n = length(bmp)
    # Initialization of the maps
    g = Dict{BitVector, UInt64}()
    h = Dict{BitVector, UInt64}()
    prev = Dict{BitVector, UInt64}()
    # The priority queue
    pq = CustomHeap()
    zero_state = falses(n)
    g[zero_state] = 0
    h[zero_state] = 0
    prev[zero_state] = 0
    update!(pq, zero_state, 0)
    # Loop
    found = false
    while !found
        cost, state = pop!(pq)
        path = trace_path(state, prev)
        k = length(path)
        partial_order!(bmp, path)
        if k == n
            found = true
            break
        end
        for i=1:n
            if state[i]
                continue
            end
            newstate = copy(state)
            newstate[i] = true
            newcost = g[state] + bonddims(bmp, k+1)
            if haskey(g, newstate)
                if newcost >= g[newstate]
                    continue
                end
            end
            g[newstate] = newcost
            prev[newstate] = i
            if !haskey(h, newstate)
                src = bmp.position[i]
                for ind=src-1:-1:k+1
                    swap!(bmp, ind)
                end
                if k < n-1
                    bdim = bonddims(bmp, k+2)
                    h[newstate] = bdim + compute_heuristic(bdim, n-k-2)
                else
                    h[newstate] = 0
                end
            end
            update!(pq, newstate, newcost+h[newstate])
        end
    end
end

"""
    exact_minimize!(bmp::BMP)

Finds the optimal variable ordering for `bmp` using an algorithm based on A*
and branch-and-bound techniques.
"""
function exact_minimize!(bmp::BMP)
    n = length(bmp)
    # Initialization of the maps
    g = Dict{BitVector, UInt64}()
    h = Dict{BitVector, UInt64}()
    lastvar = Dict{BitVector, UInt64}()
    # The upper bound
    min_vol = volume(bmp)
    min_order = copy(bmp.order)
    # The priority queue
    pq = CustomHeap()
    zero_state = falses(n)
    g[zero_state] = 0
    h[zero_state] = 0
    lastvar[zero_state] = 0
    update!(pq, zero_state, 0)
    # Loop
    found = false
    while !found
        cost, state = pop!(pq)
        # Abort if the lower bound exceeds the upper bound
        if cost + 2 >= min_vol
            found = true
            partial_order!(bmp, min_order)
            break
        end
        # Reconstruct partial variable ordering & update bounds
        path = trace_path(state, lastvar)
        k = length(path)
        min_vol = reconstruct_ordering!(bmp, path, min_vol, min_order)
        # Abort if the target state has been reached
        if k == n
            found = true
            break
        end
        # Insert successor states to the queue
        for i=1:n
            if state[i]
                continue
            end
            newstate = copy(state)
            newstate[i] = true
            newcost = g[state] + bonddims(bmp, k+1)
            if haskey(g, newstate)
                if newcost >= g[newstate]
                    # Shorter path already known, skip this state
                    continue
                end
            end
            g[newstate] = newcost
            lastvar[newstate] = i
            # Compute the lower bound estimate if necessary
            if newcost + compute_heuristic(bonddims(bmp, k+1), n-k-1) + 2 > min_vol
                continue
            end
            if !haskey(h, newstate)
                src = bmp.position[i]
                for ind=src-1:-1:k+1
                    swap!(bmp, ind)
                    vol = volume(bmp)
                    if vol < min_vol
                        min_order .= bmp.order
                    end
                end
                if k < n-1
                    bdim = bonddims(bmp, k+2)
                    h[newstate] = bdim + compute_heuristic(bdim, n-k-2)
                else
                    h[newstate] = 0
                end
            end
            # Add the successor state to the queue if necessary
            lb = h[newstate]
            if newcost + lb + 2 <= min_vol
                update!(pq, newstate, newcost + lb)
            end
        end
    end
end

"""
    sift!(bmp::BMP, n_iters::Integer)

A heuristic variable ordering optimizer that works by finding the optimal
position for one variable at a time while keeping the relative ordering of the
others fixed. `n_iters` determines how many times this process is repeated. The
default value of `n_iters` is `1`.
"""
function sift!(bmp::BMP, n_iters::Integer=1)
    n = length(bmp)
    for iter=1:n_iters
        for var_i=1:n
            pos = bmp.position[var_i]
            min_vol = volume(bmp)
            min_pos = pos
            for i=Iterators.flatten((pos-1:-1:1, 1:n-1))
                swap!(bmp, i)
                vol = volume(bmp)
                if vol < min_vol
                    min_vol = vol
                    min_pos = bmp.position[var_i]
                end
            end
            for i=n-1:-1:min_pos
                swap!(bmp, i)
            end
        end
    end
end
