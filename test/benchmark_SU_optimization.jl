using Random

function lrstep_default(mats::Matrix{Vector{UInt32}})
    U = cat(mats[1,1], mats[1,2], dims=1)
    unique!(sort!(U))
    result = Matrix{Vector{UInt32}}(undef, (2,2))
    result[1,1] = [searchsortedfirst(U, x) for x in mats[1,1]]
    result[1,2] = [searchsortedfirst(U, x) for x in mats[1,2]]
    result[2,1] = [mats[2,1][x] for x in U]
    result[2,2] = [mats[2,2][x] for x in U]
    return result
end

function lrstep_optimized(mats::Matrix{Vector{UInt32}})
    U = Dict{UInt32, UInt32}()
    sizehint!(U, length(mats[2,1]))
    L0 = [get!(U, i, length(U)+1) for i in mats[1,1]]
    L1 = [get!(U, i, length(U)+1) for i in mats[1,2]]
    result = Matrix{Vector{UInt32}}(undef, (2,2))
    result[1,1] = L0
    result[1,2] = L1
    R0 = fill(UInt32(0), length(U))
    R1 = fill(UInt32(0), length(U))
    for k in keys(U)
        ind = U[k]
        R0[ind] = mats[2,1][k]
        R1[ind] = mats[2,2][k]
    end
    result[2,1] = R0
    result[2,2] = R1
    return result
end

function rlstep_default(mats::Matrix{Vector{UInt32}})
    A = [x for x in zip(mats[2,1], mats[2,2])]
    unique!(sort!(A))
    result = Matrix{Vector{UInt32}}(undef, (2,2))
    result[1,1] = [searchsortedfirst(A, (mats[2,1][i], mats[2,2][i])) for i in mats[1,1]]
    result[1,2] = [searchsortedfirst(A, (mats[2,1][i], mats[2,2][i])) for i in mats[1,2]]
    result[2,1] = fill(UInt32(0), length(A))
    result[2,2] = fill(UInt32(0), length(A))
    for (i,pair) in enumerate(A)
        result[2,1][i] = pair[1]
        result[2,2][i] = pair[2]
    end
    return result
end

function rlstep_optimized(mats::Matrix{Vector{UInt32}})
    U = Dict{Tuple{UInt32, UInt32}, UInt32}()
    sizehint!(U, length(mats[2,1]))
    for pair in zip(mats[2,1], mats[2,2])
        get!(U, pair, length(U)+1)
    end
    result = Matrix{Vector{UInt32}}(undef, (2,2))
    result[1,1] = [get(U, (mats[2,1][i], mats[2,2][i]), 0) for i in mats[1,1]]
    result[1,2] = [get(U, (mats[2,1][i], mats[2,2][i]), 0) for i in mats[1,2]]
    result[2,1] = fill(UInt32(0), length(U))
    result[2,2] = fill(UInt32(0), length(U))
    for k in keys(U)
        ind = U[k]
        result[2,1][ind] = k[1]
        result[2,2][ind] = k[2]
    end
    return result
end

function generate_sample(chi1, chi2, chi3)
    result = Matrix{Vector{UInt32}}(undef, (2,2))
    result[1,1] = rand(1:chi2, chi1)
    result[1,2] = rand(1:chi2, chi1)
    result[2,1] = rand(1:chi3, chi2)
    result[2,2] = rand(1:chi3, chi2)
    return result
end

function verify(result1, result2)
    return all(length.(result1) .== length.(result2))
end

let
    compile_sample = generate_sample(4,4,4)
    compile_result_default = lrstep_default(compile_sample)
    compile_result_optimized = lrstep_optimized(compile_sample)
    compile_sample = generate_sample(4,4,4)
    compile_result_default = rlstep_default(compile_sample)
    compile_result_optimized = rlstep_optimized(compile_sample)
    #
    M = 10000000
    #
    lr_sample = generate_sample(M, M, M)
    lr_result_def = @time lrstep_default(lr_sample)
    lr_result_opt = @time lrstep_optimized(lr_sample)
    @show verify(lr_result_def, lr_result_opt)
    #
    rl_sample = generate_sample(M, M, M)
    rl_result_def = @time rlstep_default(rl_sample)
    rl_result_opt = @time rlstep_optimized(rl_sample)
    @show verify(rl_result_def, rl_result_opt)
end
