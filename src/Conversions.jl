include("BMP.jl")
include("BDD.jl")

function BMP(bdd::BDD)
    n = length(bdd)
    V = length(bdd.nodes)
    level_arr = Matrix{Vector{RSMInt}}(undef, (n,2))
    funcs = bdd.outputs # Outputs of the truncated BMP
    for i=1:n
        level_ids = Dict{UInt32, RSMInt}()
        level_arr[i,1] = fill(RSMInt(0), length(funcs))
        level_arr[i,2] = fill(RSMInt(0), length(funcs))
        for (j,fid) in enumerate(funcs)
            var_idx, lchild, rchild = bdd.nodes[fid]
            if var_idx != i
                # Pass-through function, convert to a redundant node
                lchild = fid
                rchild = fid
            end
            level_arr[i,1][j] = get!(level_ids, lchild, length(level_ids)+1)
            level_arr[i,2][j] = get!(level_ids, rchild, length(level_ids)+1)
        end
        funcs = fill(UInt32(0), length(level_ids))
        for (k,v) in pairs(level_ids)
            funcs[v] = k
        end
    end
    R = [bdd.nodes[fid][3] for fid in funcs]
    mats = Matrix{RowSwitchMatrix}(undef, (n,2))
    for i=1:n, j=1:2
        ncols = length(R)
        if i < n
            ncols = length(level_arr[i,j])
        end
        mats[i,j] = RowSwitchMatrix(level_arr[i,j], ncols)
    end
    return BMP(mats, R, copy(bdd.order))
end

function BDD(bmp::BMP)
end
