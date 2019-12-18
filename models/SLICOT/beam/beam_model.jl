# ===========================================================
# Beam
#
# system type: LTI system
# state dimension: 348
# input dimension: 1
# ===========================================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function beam_model()
    file = matopen(@relpath "beam.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = read(file, "B")

    # input domain
    U = BallInf([0.5], 0.3)

    # continuous LTI system
    S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U)

    return S
end
