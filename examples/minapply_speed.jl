using BinaryMatrixProducts
using Random

function compare_apply(n::Integer)
    f = generate_bmp(n, 1, rand(0:1, 2^n))
    g = generate_bmp(n, 1, rand(0:1, 2^n))
    h = rand(0:1, 4)
    println("n = $(n), k = 2")
    bmp1 = @time apply(f, g, h)
    bmp2 = @time minapply(f, g, h)
    tests_in = bitrand(n, 1000)
    @show all(evalfunc(bmp1, tests_in) .== evalfunc(bmp2, tests_in))
end

function compare_apply(n::Integer, k::Integer)
    fs = [generate_bmp(n, 1, rand(0:1, 2^n)) for _ in 1:k]
    h = rand(0:1, 2^k)
    println("n = $(n), k = $(k)")
    bmp1 = @time apply(fs, h)
    bmp2 = @time minapply(fs, h)
    tests_in = bitrand(n, 1000)
    @show all(evalfunc(bmp1, tests_in) .== evalfunc(bmp2, tests_in))
end

function build_adder(n::Integer, apply_func::Function, good_order::Bool)
    g_bit = [0, 1, 1, 0, 1, 0, 0, 1]
    g_carry = [0, 0, 0, 1, 0, 1, 1, 1]
    c = BMP(0, 2*n)
    outputs = Vector{BMP}(undef, n+1)
    for i=1:n
        xi = i
        yi = i+n
        if good_order
            xi = 2*n - 2*i + 1
            yi = 2*n - 2*i + 2
        end
        x = projbmp(xi, 2*n)
        y = projbmp(yi, 2*n)
        outputs[i] = apply_func([x, y, c], g_bit)
        c = apply_func([x, y, c], g_carry)
    end
    outputs[n+1] = c
    return joinfuncs(outputs)
end

let
    println("FORCE COMPILATION")
    compare_apply(4)
    compare_apply(4, 3)
    println("\nBINARY OPERATIONS")
    for n=4:2:16
        compare_apply(n)
    end
    println("\nTERNARY OPERATIONS")
    for n=4:2:12
        compare_apply(n, 3)
    end
    #
    build_adder(4, apply, true)
    build_adder(4, minapply, true)
    build_adder(4, apply, false)
    build_adder(4, minapply, false)
    println("\nADDERS, GOOD VARIABLE ORDERING")
    for n=4:2:16
        println("n = $(n)")
        bmp1 = @time build_adder(n, apply, true)
        bmp2 = @time build_adder(n, minapply, true)
        tests_in = bitrand(2*n, 1000)
        @show all(evalfunc(bmp1, tests_in) .== evalfunc(bmp2, tests_in))
    end
    println("\nADDERS, BAD VARIABLE ORDERING")
    for n=4:2:16
        println("n = $(n)")
        bmp1 = @time build_adder(n, apply, false)
        bmp2 = @time build_adder(n, minapply, false)
        tests_in = bitrand(2*n, 1000)
        @show all(evalfunc(bmp1, tests_in) .== evalfunc(bmp2, tests_in))
    end
end
