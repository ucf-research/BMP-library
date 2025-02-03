include("../src/BMP.jl")

using Random

function random_bmp(n::Integer)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        d = 2^(i-1)
        mats[i,1] = RowSwitchMatrix(collect(1:d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(d+1:2*d), 2*d)
    end
    return BMP_clean1(BMP(mats, rand(0:1, 2^n), collect(1:n)))
end

function all_inputs(n::Integer)
    return [val >> i & 1 for i=n-1:-1:0, val=0:2^n-1]
end

let
    n = 16
    f = random_bmp(n)
    @time BMP_save("test_bmp.txt", f)
    ff = @time BMP_load("test_bmp.txt")
    inp = all_inputs(n)
    @show all(eval(f, inp) .== eval(f, inp))
end
