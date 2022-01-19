module Moperator

using LinearAlgebra
using SparseArrays

export moperator

# This is a tensor of the operator applying on the k-th spin.
function moperator(v::Matrix, k::Int, N::Int)
    m = 1
    for j ∈ 1:N
        if j == k
            tv =v
        else
            tv = Matrix(I, 2, 2)
        end
        m = sparse(kron(m, tv))
    end
    return m
end

end
