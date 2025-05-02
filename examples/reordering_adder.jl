using BinaryMatrixProducts
using BinaryMatrixProducts: basic_exact_minimize!

function build_adder(n::Integer)
    c = BMP(0, 2*n)
    outputs = Vector{BMP}()
    g_digit = [0, 1, 1, 0, 1, 0, 0, 1]
    g_carry = [0, 0, 0, 1, 0, 1, 1, 1]
    for i=1:n
        xi = projbmp(i, 2*n)
        yi = projbmp(i+n, 2*n)
        push!(outputs, apply([xi, yi, c], g_digit))
        c = apply([xi, yi, c], g_carry)
    end
    return joinfuncs([outputs; c])
end

function time_methods(n::Integer)
    adder1 = build_adder(n)
    adder2 = BMP(copy(adder1.M), copy(adder1.R), copy(adder1.order))
    adder3 = BMP(copy(adder1.M), copy(adder1.R), copy(adder1.order))
    println("n = $(n)")
    @time exact_minimize!(adder1)
    @show volume(adder1)
    @time basic_exact_minimize!(adder2)
    @show volume(adder2)
    @time sift!(adder3)
    @show volume(adder3)
    println()
end

let
    # Note that this takes several minutes to run
    for n=8:2:12
        time_methods(n)
    end
end
