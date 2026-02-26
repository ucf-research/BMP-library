using BinaryMatrixProducts

let
    # Construct the function f(x) = x1 x2 + x3 x4 + x5 x6 + x7 x8 + x9 x10
    n = 10
    f = BMP(0, n)
    for i=1:2:n
        xl = projbmp(i, n)
        xr = projbmp(i+1, n)
        temp = apply([0, 0, 0, 1, 1, 1, 1, 1], [f, xl, xr]) # f = f + xl xr
        f = temp
    end

    # This function has a small BMP
    @show bonddims(f)
    @show volume(f), volume(f) == 5 * div(n,2)

    # The size of the BMP is sensitive to the variable ordering
    order = vcat(collect(1:2:n), collect(2:2:n))
    @show order
    reorder!(f, order)
    @show bonddims(f)
    @show volume(f), volume(f) == 3 * 2 ^ div(n,2) + div(n,2) - 2

    # We can use exact minimization to find the optimal variable ordering
    exact_minimize!(f)
    @show Vector{Int}(f.order)
    @show bonddims(f)
    @show volume(f), volume(f) == 5 * div(n,2)

    # A heuristic method called sifting can also be used, but it is not
    # guaranteed to find the optimal variable ordering
    reorder!(f, order)
    sift!(f)
    @show Vector{Int}(f.order)
    @show bonddims(f)
    @show volume(f)
end
