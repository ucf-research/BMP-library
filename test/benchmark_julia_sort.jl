using Random

function test_times(k::Integer, N::Integer)
    a = rand(1:N, (N, k))
    b = [Tuple(a[i,:]) for i=1:size(a,1)]
    println(sizeof(a), ", ", sizeof(b))
    @time a_ = sortslices(a, dims=1)
    @time b_ = sort(b)
    @show all(b_[i] == Tuple(a_[i,:]) for i=1:size(a,1))
end

let
    test_times(2, 10)
    println("k = 2, N = 1,000,000")
    test_times(2, 1000000)
    println("k = 2, N = 5,000,000")
    test_times(2, 5000000)
    println("k = 2, N = 10,000,000")
    test_times(2, 10000000)
    println("k = 3, N = 1,000,000")
    test_times(3, 1000000)
    println("k = 3, N = 10,000,000")
    test_times(3, 10000000)
    println("k = 5, N = 1,000,000")
    test_times(5, 1000000)
    println("k = 5, N = 10,000,000")
    test_times(5, 10000000)
    println("k = 10, N = 1,000,000")
    test_times(10, 1000000)
    println("k = 10, N = 10,000,000")
    test_times(10, 10000000)
end
