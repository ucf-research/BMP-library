using BinaryMatrixProducts
using BinaryMatrixProducts: apply_noclean
using BinaryMatrixProducts: minapply_noclean
using Test
using Random

@testset "APPLY methods comparison" begin
    n = 6
    all_inputs = BitMatrix(undef, (n, 2^n))
    for val=0:2^n-1
        for i=1:n
            all_inputs[i, val+1] = val >> (n-i) & 1
        end
    end
    for _ in 1:5
        bmp1 = generate_bmp(n, 1, rand(0:1, 2^n))
        bmp2 = generate_bmp(n, 1, rand(0:1, 2^n))
        htab = rand(0:1, 4)
        large1 = apply_noclean(bmp1, bmp2, htab)
        large2 = minapply_noclean(bmp1, bmp2, htab)
        @test all(bonddims(large1) .>= bonddims(large2))
        @test all(evalfunc(large1, all_inputs) .== evalfunc(large2, all_inputs))
        small1 = apply(bmp1, bmp2, htab)
        small2 = minapply(bmp1, bmp2, htab)
        @test all(bonddims(small1) .== bonddims(small2))
        @test all(evalfunc(small1, all_inputs) .== evalfunc(small2, all_inputs))
        @test all(evalfunc(small1, all_inputs) .== evalfunc(large1, all_inputs))
    end
    for _ in 1:5
        input_bmps = [
            generate_bmp(n, 1, rand(0:1, 2^n))
            generate_bmp(n, 1, rand(0:1, 2^n))
            generate_bmp(n, 1, rand(0:1, 2^n))
        ]
        htab = rand(0:1, 8)
        large1 = apply_noclean(input_bmps, htab)
        large2 = minapply_noclean(input_bmps, htab)
        @test all(bonddims(large1) .>= bonddims(large2))
        @test all(evalfunc(large1, all_inputs) .== evalfunc(large2, all_inputs))
        small1 = apply(input_bmps, htab)
        small2 = minapply(input_bmps, htab)
        @test all(bonddims(small1) .== bonddims(small2))
        @test all(evalfunc(small1, all_inputs) .== evalfunc(small2, all_inputs))
        @test all(evalfunc(small1, all_inputs) .== evalfunc(large1, all_inputs))
    end
end
