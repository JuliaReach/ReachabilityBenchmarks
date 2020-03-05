# ==================================
# Heat
#
# system type: continuous LTI system
# state dimension: 200
# input dimension: 1
# ==================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT, SparseArrays

function heat_model()
    file = matopen(@relpath "heat.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = sparse([67], [1], [1.0], size(A, 1), 1)

    # state domain
    X = Universe(200)

    # input domain
    U = BallInf([0.0], 0.5)

    # continuous LTI system
    S = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    return S
end
