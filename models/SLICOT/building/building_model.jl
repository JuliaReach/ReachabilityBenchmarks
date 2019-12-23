# ===========================================================
# Building
#
# system type: LTI system
# state dimension: 48
# input dimension: 1
# ===========================================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function building_model()
    file = matopen(@relpath "building.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = read(file, "B")

    # input domain
    U = BallInf([0.9], 0.1)

    # continuous LTI system
    S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U)

    return S
end
