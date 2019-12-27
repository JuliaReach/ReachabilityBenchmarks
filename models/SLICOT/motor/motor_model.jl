# ========================
# Motor
#
# system type: LTI system
# state dimension: 8
# input dimension: 2
# ========================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, SparseArrays

function motor_model()
    # system matrix
    I = [1, 2, 2, 3, 3, 3, 3, 4, 5, 6, 6, 7, 7, 7, 7, 8]
    J = [2, 3, 2, 1, 2, 3, 4, 1, 6, 7, 6, 5, 6, 7, 8, 5]
    V = [1, 8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0, 1.0,
         8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0]
    A = sparse(I, J, V)

    # input matrix
    B = sparse([4, 8], [1, 2], [-1.0, -1.0])

    # input domain
    U = Hyperrectangle([0.23, 0.3], [0.07, 0.1])

    # continuous LTI system
    S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U)

    return S
end
