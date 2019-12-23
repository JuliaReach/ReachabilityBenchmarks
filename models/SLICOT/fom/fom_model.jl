# ===========================================================
# FOM
#
# system type: LTI system
# state dimension: 1006
# input dimension: 1
# ===========================================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function fom_model()
    file = matopen(@relpath "fom.mat")

    # system matrix
    A = float(read(file, "A"))  # the matrix has Int entries

    # input matrix
    B = read(file, "B")

    # input domain
    U = BallInf([0.0], 1.0)

    # continuous LTI system
    S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U)

    return S
end
