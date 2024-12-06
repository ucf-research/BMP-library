include("../src/BMP_ordering.jl")

using Random

let
    heap = CustomHeap()
    for i=1:20
        label = bitrand(4)
        cost = rand(0:1000)
        println("Label: ", bitstring(label), ", Cost: ", cost)
        update!(heap, label, cost)
    end
    println("\nElements in order: ")
    while length(heap.costs) > 0
        c, s = pop!(heap)
        println("Label: ", bitstring(s), ", Cost: ", c)
    end
end
