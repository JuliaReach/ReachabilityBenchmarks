# ==================================
# PDE
#
# system type: continuous LTI system
# state dimension: 84
# input dimension: 1
# ==================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function pde_model()
    file = matopen(@relpath "pde.mat")

    # system matrix
    A = float(read(file, "A"))  # the matrix has Int entries

    # input matrix
    B = read(file, "B")

    # state domain
    X = Universe(84)

    # input domain
    U = BallInf([0.75], 0.25)

    # continuous LTI system
    S = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    return S
end
