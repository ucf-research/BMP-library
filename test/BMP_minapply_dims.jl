include("../src/BMP.jl")

using Random

function generate_function_bmp(n::Integer, f::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        d = 2^(i-1)
        mats[i,1] = RowSwitchMatrix(collect(1:d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(d+1:2*d), 2*d)
    end
    return BMP_clean1(BMP(mats, f, collect(1:n)))
end

let
    n = 16
    fvec = rand(0:1, 2^n)
    bmp_f = generate_function_bmp(n, fvec)
    gvec = rand(0:1, 2^n)
    bmp_g = generate_function_bmp(n, gvec)
    #
    mats_0 = @time BMP_apply_noclean([bmp_f.M, bmp_g.M])
    mats_1, U = @time BMP_minapply([bmp_f.M, bmp_g.M])
    @show BMP_dims(mats_0)
    @show BMP_dims(mats_1)
    @show BMP_dims(bmp_f)
    @show BMP_dims(bmp_g)
    @show BMP_dims(bmp_f) .+ BMP_dims(bmp_g)
    R0 = BMP_apply_R([bmp_f.R, bmp_g.R], [0,1,0,1])
    bmp_0 = BMP_clean1(BMP(mats_0, R0, collect(1:n)))
    R1 = BMP_minapply_R(U, [bmp_f.R, bmp_g.R], [0,1,0,1])
    bmp_1 = BMP(BMP_clean1_rl(mats_1, R1), [0,1], collect(1:n))
    @show BMP_dims(bmp_0)
    @show BMP_dims(bmp_1)
    tests = bitrand(n, 1000)
    @show all(eval(bmp_0, tests) .== eval(bmp_1, tests))
end
