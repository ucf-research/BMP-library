# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

function logged_exact_minimize!(bmp::BMP, output::IO)
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
    #
    n_expansions = 0
    n_repeats = 0
    prev_expansions = Set{BitVector}()
    # Loop
    found = false
    while !found
        cost, state = pop!(pq)
        #
        n_expansions += 1
        is_repeat = state in prev_expansions
        if is_repeat
            n_repeats += 1
        end
        push!(prev_expansions, state)
        println(
            output,
            "$(n_expansions) ",
            "$(join(string.(1 * state))) ",
            "$(count(state)) ",
            "$(is_repeat) ",
            "$(g[state]) ",
            "$(h[state]) ",
            "$(cost) ",
            "$(g[state] + h[state]) ",
            "$(min_vol) "
        )
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
                println(output, "  $(join(string.(1 * newstate))) $(newcost) $(lb) $(newcost+lb) $(min_vol)")
                update!(pq, newstate, newcost + lb)
            end
        end
    end
    println()
    println(output, "$(n_expansions) $(n_repeats) $(length(prev_expansions)) $(volume(bmp))")
end

function logged_sift!(bmp::BMP, output::IO, n_iters::Integer=1)
    n = length(bmp)
    for iter=1:n_iters
        println("Iteration: ", iter)
        for var_i=1:n
            pos = bmp.position[var_i]
            min_vol = volume(bmp)
            min_pos = pos
            println(
                output,
                "  Variable: ",
                var_i,
                ", initial position: ",
                pos,
                ", initial volume: ",
                min_vol,
                ", initial order: ",
                Vector{Int64}(bmp.order)
            )
            for i=Iterators.flatten((pos-1:-1:1, 1:n-1))
                swap!(bmp, i)
                vol = volume(bmp)
                println(
                    output,
                    "    Current position: ",
                    bmp.position[var_i],
                    ", current volume: ",
                    vol,
                    ", current order: ",
                    Vector{Int64}(bmp.order)
                )
                if vol < min_vol
                    min_vol = vol
                    min_pos = bmp.position[var_i]
                end
            end
            println(output, "  Minimum position: ", min_pos, ", minimum volume: ", min_vol)
            for i=n-1:-1:min_pos
                swap!(bmp, i)
                vol = volume(bmp)
                println(
                    output,
                    "    Current position: ",
                    bmp.position[var_i],
                    ", current volume: ",
                    vol,
                    ", current order: ",
                    Vector{Int64}(bmp.order)
                )
            end
        end
    end
end
