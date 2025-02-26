include("../src/Conversions.jl")

let
    bdd = BDD([
        BDDNode(1, 2, 3),
        BDDNode(2, 4, 5),
        BDDNode(2, 5, 6),
        BDDNode(4, 0, 0),
        BDDNode(3, 8, 7),
        BDDNode(3, 8, 9),
        BDDNode(4, 0, 1),
        BDDNode(4, 0, 0),
        BDDNode(4, 0, 1)
    ], [1], [1,2,3])
    for (i, nd) in pairs(bdd.nodes)
        tp = (nd.var_idx, nd.lchild, nd.hchild)
        println("\t$i\t$tp\t$(bdd.ref_count[i])")
    end
    println()
    bmp = BMP(bdd)
    println("Corresponding BMP size: ", BMP_volume(bmp))
    println()
    newbdd = BDD(copy(bdd.nodes), copy(bdd.outputs), copy(bdd.order))
    BDD_reduce!(newbdd)
    for (i, nd) in pairs(newbdd.nodes)
        tp = (nd.var_idx, nd.lchild, nd.hchild)
        println("\t$i\t$tp\t$(newbdd.ref_count[i])")
    end
    println()
    newbmp = BMP(newbdd)
    println("Corresponding BMP size: ", BMP_volume(newbmp))
    println()
    bdd_ = BDD(newbmp)
    for (i, nd) in pairs(bdd_.nodes)
        tp = (nd.var_idx, nd.lchild, nd.hchild)
        println("\t$i\t$tp\t$(bdd_.ref_count[i])")
    end
    println()
    #
    n = 3
    test_input = BitArray(undef, (n, 2^n))
    for j=1:n, i=0:2^n-1
        test_input[j, i+1] = i >> (n-j) & 1
    end
    @show all(eval(bdd, test_input) .== eval(newbdd, test_input))
    @show all(eval(bdd, test_input) .== eval(bmp, test_input))
    @show all(eval(bdd, test_input) .== eval(newbmp, test_input))
    @show all(eval(bdd, test_input) .== eval(bdd_, test_input))
end
