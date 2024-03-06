# ==================================
# FOM
#
# system type: continuous LTI system
# state dimension: 1006
# input dimension: 1
# ==================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function fom_model()
    file = matopen(@current_path "fom.mat")

    # system matrix
    A = float(read(file, "A"))  # the matrix has Int entries

    # input matrix
    B = read(file, "B")

    # state domain
    X = Universe(1006)

    # input domain
    U = BallInf([0.0], 1.0)

    # continuous LTI system
    S = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    return S
end
