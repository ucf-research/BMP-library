include("../src/RowSwitchMatrix.jl")

let
    a = RowSwitchMatrix([1,2,3], 3)
    b = RowSwitchMatrix([2,1], 2)
    display(RSM_matrix(a))
    display(RSM_matrix(b))
    display(RSM_matrix(RSM_kron([a,b])))
    display(RSM_matrix(RSM_join([a,b])))
end
