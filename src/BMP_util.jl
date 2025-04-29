function generate_bmp(n::Integer, m::Integer, f::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    d = m
    for i=1:n
        mats[i,1] = RowSwitchMatrix(collect(1:2:2*d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(2:2:2*d), 2*d)
        d *= 2
    end
    return clean1(BMP(mats, f, collect(1:n)))
end

