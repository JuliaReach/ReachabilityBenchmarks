# ==================================
# Building
#
# system type: continuous LTI system
# state dimension: 48
# input dimension: 1
# ==================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function building_model()
    file = matopen(@current_path "building.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = read(file, "B")

    # state domain
    X = Universe(48)

    # input domain
    U = BallInf([0.9], 0.1)

    # continuous LTI system
    S = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    return S
end
