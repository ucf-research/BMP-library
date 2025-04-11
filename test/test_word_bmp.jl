include("../src/WordBMP.jl")

using Random

function generate_function_bmp(n::Integer, m::Integer, f::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        d = 2^(i-1) * m
        mats[i,1] = RowSwitchMatrix(collect(1:2:2*d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(2:2:2*d), 2*d)
    end
    return BMP_clean1(BMP(mats, f, collect(1:n)))
end

function verify_conversions(n::Integer, m::Integer, ki::Integer, ko::Integer)
    f = rand(0:1, m * 2^n)
    bmp1 = @time generate_function_bmp(n, m, f)
    wbmp = @time WordBMP(bmp1, ki, ko)
    bmp2 = @time BMP(wbmp)
    tests_in = bitrand(n, 1000)
    # @show all(eval(bmp1, tests_in) .== f)
    @show all(eval(bmp1, tests_in) .== eval(wbmp, tests_in))
    @show all(eval(bmp1, tests_in) .== eval(bmp2, tests_in))
    @show BMP_volume(bmp1), BMP_volume(bmp2)
end

let
    verify_conversions(12, 8, 3, 4)
    verify_conversions(12, 8, 3, 4)
end
