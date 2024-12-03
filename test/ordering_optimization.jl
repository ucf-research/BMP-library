include("../src/BMP_ordering.jl")

using Random
using Combinatorics

function generate_function_bmp(n::Integer, f::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        d = 2^(i-1)
        mats[i,1] = RowSwitchMatrix(collect(1:d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(d+1:2*d), 2*d)
    end
    return BMP_clean1(BMP(mats, f, collect(1:n)))
end

function brute_force!(bmp::BMP)
    n = length(bmp)
    min_vol = BMP_volume(bmp)
    min_order = copy(bmp.order)
    for perm in permutations(1:n)
        partial_order!(bmp, perm)
        vol = BMP_volume(bmp)
        if vol < min_vol
            min_vol = vol
            min_order = copy(perm)
        end
    end
    partial_order!(bmp, min_order)
    return (min_vol, min_order)
end

function Astar!(bmp::BMP)
    basic_exact_minimize!(bmp)
    return (BMP_volume(bmp), copy(bmp.order))
end

function Astar_BB!(bmp::BMP)
    exact_minimize!(bmp)
    return (BMP_volume(bmp), copy(bmp.order))
end

function dynamic!(bmp::BMP)
    sift!(bmp)
    return (BMP_volume(bmp), copy(bmp.order))
end

function benchmark_all(n::Integer)
    f = rand(0:1, 2^n)
    bmp1 = generate_function_bmp(n, f)
    bmp2 = generate_function_bmp(n, f)
    bmp3 = generate_function_bmp(n, f)
    bmp4 = generate_function_bmp(n, f)
    v1, ord1 = @time brute_force!(bmp1)
    v2, ord2 = @time Astar!(bmp2)
    v3, ord3 = @time Astar_BB!(bmp3)
    v4, ord4 = @time dynamic!(bmp4)
    println("n = ", n)
    @show v1, Vector{Int64}(ord1)
    @show v2, Vector{Int64}(ord2)
    @show v3, Vector{Int64}(ord3)
    @show v4, Vector{Int64}(ord4)
    println()
end

function benchmark_smart(n::Integer)
    f = rand(0:1, 2^n)
    bmp1 = generate_function_bmp(n, f)
    bmp2 = generate_function_bmp(n, f)
    v1, ord1 = @time Astar!(bmp1)
    v2, ord2 = @time Astar_BB!(bmp2)
    v4, ord4 = @time dynamic!(bmp4)
    println("n = ", n)
    @show v1, Vector{Int64}(ord1)
    @show v2, Vector{Int64}(ord2)
    @show v4, Vector{Int64}(ord4)
    println()
end

let
    benchmark_all(4)
    benchmark_all(5)
    benchmark_all(6)
    benchmark_all(7)
    benchmark_all(8)
    benchmark_all(9)
    benchmark_all(10)
end
