include("../src/BMP.jl")

function build_bmp(order)
    n = length(order)
    bmp = BMP(0, order)
    for i=1:2:n
        x1 = BMP_bitline(i, order)
        x2 = BMP_bitline(i+1, order)
        t = BMP_apply(x1, x2, [0,0,0,1])
        temp = BMP_apply(bmp, t, [0,1,1,1])
        bmp = temp
    end
    return bmp
end

function compare_orders(k::Integer)
    n = 2*k
    ord1 = collect(1:n)
    ord2 = collect(Iterators.flatten([1:2:n, 2:2:n]))
    bmp1 = build_bmp(ord1)
    bmp2 = build_bmp(ord2)
    vol1_th = 5 * k
    vol2_th = 3 * 2^k + k - 2
    @show vol1_th, vol2_th
    @show BMP_volume(bmp1) == vol1_th
    @show BMP_volume(bmp2) == vol2_th
end

let
    compare_orders(4)
    compare_orders(5)
    compare_orders(6)
    compare_orders(7)
    compare_orders(8)
    compare_orders(9)
    compare_orders(10)
end
