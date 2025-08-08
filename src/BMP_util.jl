# This file is a part of BMP-library. License is Apache 2.0: https://julialang.org/license

"""
    generate_bmp(n::Integer, m::Integer, f::AbstractArray)

Generates a BMP of `n` inputs and `m` outputs directly from its truth table given
in `f`. This function has necessarily exponential complexity and it should only
be used for small systems for testing purposes.
"""
function generate_bmp(n::Integer, m::Integer, f::AbstractArray)
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    d = m
    for i=1:n
        mats[i,1] = RowSwitchMatrix(collect(1:2:2*d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(2:2:2*d), 2*d)
        d *= 2
    end
    return clean1(BMP(mats, f, collect(1:n)))
end

"""
    check_equivalence(bmp1::BMP, bmp2::BMP)

Returns `true` if the arguments `bmp1` and `bmp2` represent the same function,
`false` otherwise. This function uses `minapply` to compute the XOR of the two
BMPs and checks if the result is identically zero.
"""
function check_equivalence(bmp1::BMP, bmp2::BMP)
    bmp = minapply_noclean([0, 1, 1, 0], bmp1, bmp2)
    return all(bmp.R .== 0)
end

