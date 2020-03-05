# ==================================
# CD Player
#
# system type: continuous LTI system
# state dimension: 120
# input dimension: 2
# ==================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function cdplayer_model()
    file = matopen(@relpath "cdplayer.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = read(file, "B")

    # state domain
    X = Universe(120)

    # input domain
    U = BallInf([0.0, 0.0], 1.0)

    # continuous LTI system
    S = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    return S
end
