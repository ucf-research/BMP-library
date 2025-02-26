include("../src/BMP_ordering.jl")

function generate_function_bmp(n::Integer, f::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        d = 2^(i-1)
        mats[i,1] = RowSwitchMatrix(collect(1:d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(d+1:2*d), 2*d)
    end
    return BMP_clean1(BMP(mats, f, collect(1:n)))
end

function Astar!(bmp::BMP)
    basic_exact_minimize!(bmp)
    return (BMP_volume(bmp), copy(bmp.order))
end

function Astar_BB!(bmp::BMP)
    exact_minimize!(bmp)
    return (BMP_volume(bmp), copy(bmp.order))
end

function dynamic!(bmp::BMP)
    sift!(bmp)
    return (BMP_volume(bmp), copy(bmp.order))
end

function benchmark_smart(n::Integer)
    f = rand(0:1, 2^n)
    bmp1 = generate_function_bmp(n, f)
    bmp2 = generate_function_bmp(n, f)
    bmp4 = generate_function_bmp(n, f)
    v1, ord1 = Astar!(bmp1)
    v2, ord2 = Astar_BB!(bmp2)
    v4, ord4 = dynamic!(bmp4)
end

function build_adder(n::Integer)
    outputs = BMP[]
    f_digit = [0,1,1,0,1,0,0,1]
    f_carry = [0,0,0,1,0,1,1,1]
    c = BMP(0, 2*n)
    for i=1:n
        xi = BMP_bitline(i, 2*n)
        yi = BMP_bitline(i+n, 2*n)
        push!(outputs, BMP_minapply([xi, yi, c], f_digit))
        c_ = BMP_minapply([xi, yi, c], f_carry)
        c = c_
    end
    # return BMP_join([outputs; c])
    return BMP_join(outputs)
end

function eval_adder(n, x, y, adder)
    inp = BitVector(undef, 2*n)
    inp[1:n] .= [x >> i & 1 for i=0:n-1]
    inp[n+1:2*n] .= [y >> i & 1 for i=0:n-1]
    out = eval(adder, inp)
    out_int = sum(out[i+1] * 2^i for i=0:n)
    result = x + y == out_int
    @show result, x, y, out_int
    return result
end

let
    # compile functions
    benchmark_smart(4)
    # test
    n = parse(Int, ARGS[1])
    adder1 = build_adder(n)
    adder2 = build_adder(n)
    adder3 = build_adder(n)
    @show BMP_volume(adder1)
    v1, ord1 = @time Astar!(adder1)
    v2, ord2 = @time Astar_BB!(adder2)
    v3, ord3 = @time dynamic!(adder3)
    v4, ord4 = @time Astar_BB!(adder3)
    @show v1, Vector{Int64}(ord1)
    @show v2, Vector{Int64}(ord2)
    @show v3, Vector{Int64}(ord3)
    @show v4, Vector{Int64}(ord4)
end
