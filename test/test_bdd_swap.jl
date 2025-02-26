include("../src/Conversions.jl")

using Random

function generate_function_bmp(n::Integer, f::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        d = 2^(i-1)
        mats[i,1] = RowSwitchMatrix(collect(1:d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(d+1:2*d), 2*d)
    end
    return BMP_clean1(BMP(mats, f, collect(1:n)))
end

function generate_bdd(n::Integer)
    bmp = generate_function_bmp(n, rand(0:1, 2^n))
    bdd = BDD(bmp)
    return bdd
end

function print_bdd(bdd::BDD)
    for node_ref in axes(bdd.nodes,1)
        nd = bdd.nodes[node_ref]
        tp_ = (nd.var_idx, nd.lchild, nd.hchild)
        tp = Int64.(tp_)
        print("\t$node_ref\t$tp\t$(bdd.ref_count[node_ref])\t")
        if node_ref in bdd.levels[nd.var_idx]
            print("\t")
        else
            print("X\t")
        end
        if node_ref in bdd.gc_refs
            print("X\n")
        else
            print("\n")
        end
    end
end

function random_swaps(n::Integer, n_swaps::Integer)
    bdd = generate_bdd(n)
    flag = true
    BDD_reduce!(bdd)
    swaps = rand(1:n-1, n_swaps)
    tests_in = bitrand((n, 1000))
    tests_out = eval(bdd, tests_in)
    @show bdd.order
    @show BDD_volume(bdd)
    # print_bdd(bdd)
    println()
    for swap_i in swaps
        BDD_swap!(bdd, swap_i)
        @show bdd.order
        @show vol1 = BDD_volume(bdd)
        @show match1 = all(tests_out .== eval(bdd, tests_in))
        # print_bdd(bdd)
        BDD_reduce!(bdd)
        @show vol2 = BDD_volume(bdd)
        @show match2 = all(tests_out .== eval(bdd, tests_in))
        println()
        flag = flag && vol1 == vol2 && match1 && match2
    end
    @show flag
end

let
    random_swaps(10, 30)
end
