using BinaryMatrixProducts
using Test
using Random

function compose_table(
    n::Integer,
    f::Vector{<:Integer},
    var::Integer,
    g::Vector{<:Integer}
)
    result = Vector{Int64}(undef, 2^n)
    for val=0:2^n-1
        b = g[val+1]
        mask = (2^n - 1) ⊻ (1 << (n-var))
        fval = (val & mask) | (b << (n-var))
        result[val+1] = f[fval+1]
    end
    return result
end

@testset "BMP composition" begin
    n = 8
    f = rand(0:1, 2^n)
    g = rand(0:1, 2^n)
    bmp_f = generate_bmp(n, 1, f)
    bmp_g = generate_bmp(n, 1, g)
    for vi=1:n
        h = compose_table(n, f, vi, g)
        bmp0 = generate_bmp(n, 1, h)
        bmp1 = compose(bmp_f, vi, bmp_g)
        tests_in = bitrand(n, 1000)
        @test all(evalfunc(bmp0, tests_in) .== evalfunc(bmp1, tests_in))
    end
end
