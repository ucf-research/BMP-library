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
    exact_minimize!(bmp)
    return (BMP_volume(bmp), copy(bmp.order))
end

let
    n = 5
    f = rand(0:1, 2^n)
    bmp1 = generate_function_bmp(n, f)
    bmp2 = generate_function_bmp(n, f)
    v1, ord1 = @time brute_force!(bmp1)
    v2, ord2 = @time Astar!(bmp2)
    @show v1, ord1
    @show v2, ord2
end
