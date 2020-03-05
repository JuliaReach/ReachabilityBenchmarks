# ==================================
# International Space Station (ISS)
#
# system type: continuous LTI system
# state dimension: 270
# input dimension: 3
# ==================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT

function iss_model()
    file = matopen(@relpath "iss.mat")

    # system matrix
    A = read(file, "A")

    # input matrix
    B = read(file, "B")

    # state domain
    X = Universe(270)

    # input domain
    U = Hyperrectangle([0.05, 0.9, 0.95], [0.05, 0.1, 0.05])

    # continuous LTI system
    S = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    return S
end
