include("../src/BDD.jl")

let
    bdd = BDD([
        (1, 2, 3),
        (2, 4, 5),
        (2, 5, 6),
        (4, 0, 0),
        (3, 8, 7),
        (3, 8, 9),
        (4, 0, 1),
        (4, 0, 0),
        (4, 0, 1)
    ], [1], [1,2,3])
    for (i, nd) in enumerate(bdd.nodes)
        nn = Tuple{Int64, Int64, Int64}(nd)
        println("\t$i\t$nn")
    end
    println()
    newbdd = BDD_reduce(bdd)
    for (i, nd) in enumerate(newbdd.nodes)
        nn = Tuple{Int64, Int64, Int64}(nd)
        println("\t$i\t$nn")
    end
    println()
    n = 3
    test_input = BitArray(undef, (n, 2^n))
    for j=1:n, i=0:2^n-1
        test_input[j, i+1] = i >> (n-j) & 1
    end
    @show all(eval(bdd, test_input) .== eval(newbdd, test_input))
end
