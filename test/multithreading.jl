include("../src/BMP.jl")

using Base.Threads
using Random

# Decimation eval
function deceval(bmp::BareBMP, x::BitArray, R::Vector{<:Integer}, order::Vector{<:Integer})::BitArray
    n = size(bmp, 1)
    m = length(bmp[1,1].rows)
    n_samps = div(length(x), size(x, 1))
    x_ = reshape(x, (n, n_samps))
    result = BitArray(undef, (m, n_samps))
    for j=1:n_samps
        mats = [bmp[i,x_[order[i],j]+1] for i=1:length(order)]
        while length(mats) > 1
            L = div(length(mats) + 1, 2)
            nmats = Vector{RowSwitchMatrix}(undef, L)
            Threads.@threads for k=1:L
                kr = 2*k
                kl = kr-1
                if kr > length(mats)
                    nmats[k] = mats[kl]
                else
                    nmats[k] = RSM_mult(mats[kl], mats[kr])
                end
            end
            mats = nmats
        end
        result[:,j] .= R[mats[1].rows]
    end
    shape = [i for i in size(x)]
    shape[1] = m
    return reshape(result, Tuple(shape))
end

function deceval(bmp::BMP, x::BitArray)::BitArray
    return deceval(bmp.M, x, bmp.R, bmp.order)
end

# Multiple inputs parallelized
function pareval_in(bmp::BareBMP, x::BitArray, R::Vector{<:Integer}, order::Vector{<:Integer})::BitArray
    n = size(bmp, 1)
    m = length(bmp[1,1].rows)
    n_samps = div(length(x), size(x, 1))
    x_ = reshape(x, (n, n_samps))
    result = BitArray(undef, (m, n_samps))
    BUCKET_SIZE = 128
    n_buckets = div(n_samps-1, BUCKET_SIZE) + 1
    n_threads = Threads.nthreads()
    mats = [fill(RSMInt(0), m) for _ in 1:n_threads]
    Threads.@threads for bi=0:n_buckets-1
        b_start = bi * BUCKET_SIZE + 1
        b_end = min((bi+1) * BUCKET_SIZE, n_samps)
        mat = mats[Threads.threadid()]
        for j=b_start:b_end
            mat .= bmp[1, x_[order[1],j] + 1].rows
            for i=2:n
                rs = bmp[i, x_[order[i],j] + 1].rows
                for k=1:m
                    mat[k] = rs[mat[k]]
                end
            end
            for k=1:m
                result[k,j] = R[mat[k]]
            end
        end
    end
    shape = [i for i in size(x)]
    shape[1] = m
    return reshape(result, Tuple(shape))
end

function pareval_in(bmp::BMP, x::BitArray)::BitArray
    return pareval_in(bmp.M, x, bmp.R, bmp.order)
end

# Multiple outputs parallelized
function pareval_out(bmp::BareBMP, x::BitArray, R::Vector{<:Integer}, order::Vector{<:Integer})::BitArray
    n = size(bmp, 1)
    m = length(bmp[1,1].rows)
    n_samps = div(length(x), size(x, 1))
    x_ = reshape(x, (n, n_samps))
    result = BitArray(undef, (m, n_samps))
    output_bits = [BitArray(undef, n_samps) for _ in 1:m]
    Threads.@threads for out_i=1:m
        for j=1:n_samps
            ind = out_i
            for i=1:n
                ind = bmp[i, x_[order[i],j]+1].rows[ind]
            end
            output_bits[out_i][j] = R[ind]
        end
    end
    for out_i=1:m
        result[out_i,:] .= output_bits[out_i]
    end
    shape = [i for i in size(x)]
    shape[1] = m
    return reshape(result, Tuple(shape))
end

function pareval_out(bmp::BMP, x::BitArray)::BitArray
    return pareval_out(bmp.M, x, bmp.R, bmp.order)
end

# Benchmark
function generate_function_bmp(n::Integer, f::Vector{<:Integer})
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n
        d = 2^(i-1)
        mats[i,1] = RowSwitchMatrix(collect(1:d), 2*d)
        mats[i,2] = RowSwitchMatrix(collect(d+1:2*d), 2*d)
    end
    return BMP_clean1(BMP(mats, f, collect(1:n)))
end

function generate_multifunction_bmp(n::Integer, m::Integer)
    all_bmp = [generate_function_bmp(n, rand(0:1, 2^n)) for i=1:m]
    return BMP_join(all_bmp)
end

function benchmark(n, m, n_tests)
    @show n, m, n_tests
    bmp = generate_multifunction_bmp(n, m)
    tests = bitrand(n, n_tests)
    tests_out = Vector{BitArray}(undef, 4)
    tests_out[1] = @time eval(bmp, tests)
    tests_out[2] = @time deceval(bmp, tests)
    tests_out[3] = @time pareval_in(bmp, tests)
    tests_out[4] = @time pareval_out(bmp, tests)
    for i=2:4
        @show all(tests_out[1] .== tests_out[i])
    end
    println()
end

let
    print("""
        eval: Standard evaluation
        deceval: Decimation evaluation
        pareval_in: Inputs parallelized
        pareval_out: Outputs parallelized
    \n"""
    )
    benchmark(8, 4, 100)
    benchmark(8, 4, 100)
    benchmark(16, 4, 1000)
    benchmark(16, 4, 10000)
    benchmark(16, 4, 100000)
    benchmark(16, 8, 1000)
    benchmark(16, 8, 10000)
    benchmark(16, 8, 100000)
    benchmark(16, 16, 1000)
    benchmark(16, 16, 10000)
    benchmark(16, 16, 100000)
end

#=
using Profile

function test_pareval_in(n, m, n_tests)
    bmp = generate_multifunction_bmp(n, m)
    tests = bitrand(n, n_tests)
    tests_out = pareval_in(bmp, tests)
    return tests_out
end

let
    test_pareval_in(8, 1, 100)
    @profile test_pareval_in(16, 16, 100_000)
    Profile.print()
end
=#
