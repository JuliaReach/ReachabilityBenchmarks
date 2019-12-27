# ===========================================================
# PDE
#
# system type: LTI system
# state dimension: 84
# input dimension: 1
# ===========================================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function pde_model()
    file = matopen(@relpath "pde.mat")

    # system matrix
    A = float(read(file, "A"))  # the matrix has Int entries

    # input matrix
    B = read(file, "B")

    # input domain
    U = BallInf([0.75], .25)

    # continuous LTI system
    S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U)

    return S
end
