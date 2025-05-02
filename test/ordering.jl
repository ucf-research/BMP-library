using BinaryMatrixProducts
using BinaryMatrixProducts: get_swap_position
using Test
using Random

@testset "Permutation generation" begin
    for n=4:8
        N = prod(1:n)
        a = collect(1:n)
        perms = Set{Vector{Int64}}()
        for i=1:N
            push!(perms, copy(a))
            p = get_swap_position(n, i)
            temp = a[p]
            a[p] = a[p+1]
            a[p+1] = temp
        end
        @test length(perms) == N
    end
end

@testset "BMP ordering optimization" begin
    for n=4:8
        f1 = generate_bmp(n, 1, rand(0:1, 2^n))
        init_vol = volume(f1)
        f2 = BMP(copy(f1.M), copy(f1.R), copy(f1.order))
        f3 = BMP(copy(f1.M), copy(f1.R), copy(f1.order))
        brute_force!(f1)
        exact_minimize!(f2)
        sift!(f3)
        @test volume(f1) == volume(f2)
        @test volume(f1) <= init_vol
        @test volume(f3) <= init_vol
    end
end
