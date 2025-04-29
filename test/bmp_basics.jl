using BinaryMatrixProducts
using Test

@testset "Basic BMP operations" begin
    n = 8
    # Generate BMP from truth table
    all_inputs = BitMatrix(undef, (n, 2^n))
    all_outputs = BitMatrix(undef, (1, 2^n))
    for val=0:2^n-1
        for i=1:n
            all_inputs[i, val+1] = val >> (n-i) & 1
        end
        fval = 0
        for i=2:2:n
            fval = fval | ((val >> (n-i) & 3) == 3)
        end
        all_outputs[1, val+1] = fval
    end
    f = Vector{Int64}(vec(all_outputs))
    bmp0 = generate_bmp(n, 1, f)
    @test all(evalfunc(bmp0, all_inputs) .== all_outputs)
    # Generate the same BMP using APPLY
    g_and = [0,0,0,1]
    g_or = [0,1,1,1]
    bmp1 = BMP(0, n)
    for i=1:2:n
        x1 = projbmp(i, n)
        x2 = projbmp(i+1, n)
        t = apply(x1, x2, g_and)
        bmp1 = apply(bmp1, t, g_or)
    end
    @test all(evalfunc(bmp0, all_inputs) .== evalfunc(bmp1, all_inputs))
    # Effect of the variable ordering
    vorder = [1:2:n; 2:2:n]
    bmp2 = BMP(0, vorder)
    for i=1:2:n
        x1 = projbmp(i, vorder)
        x2 = projbmp(i+1, vorder)
        t = apply(x1, x2, g_and)
        bmp2 = apply(bmp2, t, g_or)
    end
    @test all(evalfunc(bmp0, all_inputs) .== evalfunc(bmp2, all_inputs))
    n2 = div(n, 2)
    @test volume(bmp1) == 5 * n2
    @test volume(bmp2) == 3 * 2^n2 + n2 - 2
end
