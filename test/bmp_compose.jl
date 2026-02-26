using BinaryMatrixProducts
using Test
using Random

function compose_table(
    n::Integer,
    f::Vector{<:Integer},
    vars,
    subs
)
    result = Vector{Int64}(undef, 2^n)
    for val=0:2^n-1
        fval = val
        for (v, g) in zip(vars, subs)
            b = g[val+1]
            mask = ((fval >> (n-v) & 1) ⊻ b) << (n-v)
            fval = fval ⊻ mask
        end
        result[val+1] = f[fval+1]
    end
    return result
end

@testset "BMP composition" begin
    n = 8
    n_tests = 10
    for _ in 1:n_tests
        f = rand(0:1, 2^n)
        subs = ntuple(i -> rand(0:1, 2^n), 3)
        var_perm = randperm(n)
        vars = ntuple(i -> var_perm[i], 3)
        h = compose_table(n, f, vars, subs)
        #
        var_order = randperm(n)
        bmp_f = generate_bmp(n, 1, f)
        reorder!(bmp_f, var_order)
        sub_bmps = ntuple(i -> generate_bmp(n, 1, subs[i]), 3)
        for i=1:3
            reorder!(sub_bmps[i], var_order)
        end
        bmp_h = generate_bmp(n, 1, h)
        reorder!(bmp_h, var_order)
        bmp = compose(bmp_f, vars, sub_bmps)
        @test check_equivalence(bmp_h, bmp)
    end
end
