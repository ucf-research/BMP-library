using BinaryMatrixProducts

let
    # Generate the projection BMPs
    x1 = projbmp(1, 3) # First of three variables
    x2 = projbmp(2, 3)
    x3 = projbmp(3, 3)

    # Construct functions out of these using APPLY
    g = apply(x1, x3, [0, 0, 0, 1]) # g(x) = x1 & x3
    h = apply([x1, x2, x3], [0, 0, 0, 0, 0, 0, 1, 0]) # h(x) = x1 & x2 & ~x3
    f = apply([g, h], [0, 1, 1, 1]) # f(x) = g(x) | h(x)

    # Check the outputs of these functions for a few inputs
    @show evalfunc(g, [1, 0, 1]) # = 1
    @show evalfunc(h, [1, 1, 0]) # = 1
    @show evalfunc(f, [1, 1, 1]) # = 1

    # apply uses the direct product method
    # For the direct sum method, use minapply
    g_ = minapply(x1, x3, [0, 0, 0, 1])
    h_ = minapply([x1, x2, x3], [0, 0, 0, 0, 0, 0, 1, 0])
    f_ = minapply([g, h], [0, 1, 1, 1])

    # The resulting BMPs are equivalent
    # Check this from the bond dimensions
    @show bonddims(g), bonddims(g_)
    @show bonddims(h), bonddims(h_)
    @show bonddims(f), bonddims(f_)

    # Check that both versions have the same outputs
    all_inputs = BitMatrix([
        0 0 0 0 1 1 1 1;
        0 0 1 1 0 0 1 1;
        0 1 0 1 0 1 0 1
    ])
    @show evalfunc(g, all_inputs)
    @show all(evalfunc(g, all_inputs) .== evalfunc(g_, all_inputs))
    @show evalfunc(h, all_inputs)
    @show all(evalfunc(h, all_inputs) .== evalfunc(h_, all_inputs))
    @show evalfunc(f, all_inputs)
    @show all(evalfunc(f, all_inputs) .== evalfunc(f_, all_inputs))
end
