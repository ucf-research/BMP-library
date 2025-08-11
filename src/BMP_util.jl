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

function check_equivalence(
    mats1::BareBMP,
    R1::AbstractArray,
    mats2::BareBMP,
    R2::AbstractArray
)
    if size(mats1,1) != size(mats2,1)
        return false
    end
    _, U = minapply_mats((mats1, mats2))
    R = minapply_term([0, 1, 1, 0], U, (R1, R2))
    return all(v -> v == 0, R)
end

"""
    check_equivalence(bmp1::BMP, bmp2::BMP)

Returns `true` if the arguments `bmp1` and `bmp2` represent the same function,
`false` otherwise. This function uses `minapply` to compute the XOR of the two
BMPs and checks if the result is identically zero.
"""
function check_equivalence(bmp1::BMP, bmp2::BMP)
    return check_equivalence(bmp1.M, bmp1.R, bmp2.M, bmp2.R)
end
