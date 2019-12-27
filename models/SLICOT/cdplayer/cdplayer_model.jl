# ===========================================================
# CD Player
#
# system type: LTI system
# state dimension: 120
# input dimension: 2
# ===========================================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function cdplayer_model()
    file = matopen(@relpath "cdplayer.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = read(file, "B")

    # input domain
    U = BallInf([0.0, 0.0], 1.0)

    # continuous LTI system
    S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U)

    return S
end
