using BinaryMatrixProducts
using Test
using Random

function compute_tab(htab, fns)
    result = fill(0, length(fns[1]))
    for i in axes(result, 1)
        val = 0
        for f in fns
            val = 2 * val + f[i]
        end
        result[i] = htab[val+1]
    end
    return result
end

@testset "APPLY methods comparison" begin
    n = 8
    for _ in 1:5
        f1 = rand(0:1, 2^n)
        f2 = rand(0:1, 2^n)
        htab = rand(0:1, 4)
        h = compute_tab(htab, (f1, f2))
        bmp1 = generate_bmp(n, 1, f1)
        bmp2 = generate_bmp(n, 1, f2)
        bmp_h = generate_bmp(n, 1, h)
        large1 = apply(htab, bmp1, bmp2; noclean=true)
        large2 = minapply(htab, bmp1, bmp2; noclean=true)
        @test all(bonddims(large1) .>= bonddims(large2))
        @test check_equivalence(large1, bmp_h)
        @test check_equivalence(large1, large2)
        small1 = apply(htab, bmp1, bmp2)
        small2 = minapply(htab, bmp1, bmp2)
        @test all(bonddims(small1) .== bonddims(small2))
        @test check_equivalence(small1, bmp_h)
        @test check_equivalence(small1, small2)
    end
    for _ in 1:5
        f1 = rand(0:1, 2^n)
        f2 = rand(0:1, 2^n)
        f3 = rand(0:1, 2^n)
        htab = rand(0:1, 8)
        h = compute_tab(htab, (f1, f2, f3))
        input_bmps = (
            generate_bmp(n, 1, f1),
            generate_bmp(n, 1, f2),
            generate_bmp(n, 1, f3)
        )
        bmp_h = generate_bmp(n, 1, h)
        large1 = apply(htab, input_bmps; noclean=true)
        large2 = minapply(htab, input_bmps; noclean=true)
        @test all(bonddims(large1) .>= bonddims(large2))
        @test check_equivalence(large1, bmp_h)
        @test check_equivalence(large1, large2)
        small1 = apply(htab, input_bmps)
        small2 = minapply(htab, input_bmps)
        @test all(bonddims(small1) .== bonddims(small2))
        @test check_equivalence(small1, bmp_h)
        @test check_equivalence(small1, small2)
    end
end
