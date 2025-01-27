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
    ], [1,2,3])
    for (i, n) in enumerate(bdd.nodes)
        nn = Tuple{Int64, Int64, Int64}(n)
        println("\t$i\t$nn")
    end
    println()
    BDD_reduce!(bdd)
    for (i, n) in enumerate(bdd.nodes)
        nn = Tuple{Int64, Int64, Int64}(n)
        println("\t$i\t$nn")
    end
end
