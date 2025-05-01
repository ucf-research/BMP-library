using BinaryMatrixProducts
using Test
using Random

@testset "Miscellaneous BMP operations" begin
    n = 16
    f = rand(0:1, 2^n)
    bmp0 = generate_bmp(n, 1, f)
    n_tests = 1000
    tests0 = bitrand(n, n_tests)
    #
    bmp1 = insert_var(bmp0)
    tests1 = vcat(tests0, bitrand(1, n_tests))
    @test all(evalfunc(bmp0, tests0) .== evalfunc(bmp1, tests1))
    #
    for p=1:n+1
        bmp2 = insert_var(bmp0, p)
        tests2 = vcat(tests0, bitrand(1, n_tests))
        @test all(evalfunc(bmp0, tests0) .== evalfunc(bmp2, tests2))
    end
    #
    for vi=1:n, val=0:1
        bmp3 = restrict(bmp0, vi, val)
        tests0_ = copy(tests0)
        tests0_[vi,:] .= val
        @test all(evalfunc(bmp0, tests0_) .== evalfunc(bmp3, tests0))
    end
    #
    for vi=1:n, val=0:1
        bmp4 = erase_var(bmp0, vi, val)
        tests0_ = copy(tests0)
        tests0_[vi,:] .= val
        tests4 = vcat(tests0[1:vi-1,:], tests0[vi+1:n,:])
        @test all(evalfunc(bmp0, tests0_) .== evalfunc(bmp4, tests4))
    end
    #
    n_shuffles = 10
    tests_out = evalfunc(bmp0, tests0)
    for _ in 1:n_shuffles
        order = randperm(n)
        reorder!(bmp0, order)
        @test all(evalfunc(bmp0, tests0) .== tests_out)
    end
end

