include("../src/BMP.jl")

let
    f = [x3 & x2 | x1 for x3=0:1, x2=0:1, x1=0:1]
    x1 = BMP_bitline(1, 3)
    x2 = BMP_bitline(2, 3)
    x3 = BMP_bitline(3, 3)
    bmp1 = @time BMP_apply(x2, x3, [0,0,0,1])
    bmp1 = @time BMP_apply(x1, bmp1, [0,1,1,1])
    bmp_outputs = [eval(bmp1, [x1, x2, x3]) for x3=0:1, x2=0:1, x1=0:1]
    @show all(u == v[1] for (u,v) in zip(f, bmp_outputs))
end
