# ==================================
# Beam
#
# system type: continuous LTI system
# state dimension: 348
# input dimension: 1
# ==================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function beam_model()
    file = matopen(@current_path "beam.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = read(file, "B")

    # state domain
    X = Universe(348)

    # input domain
    U = BallInf([0.5], 0.3)

    # continuous LTI system
    S = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    return S
end
