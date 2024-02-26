# ==============================================================================
# Five-dimensional system
#
# system type: continuous LTI system
# state dimension: 5
# input dimension: 1
#
# This five-dimensional system is taken from Example 4.1 in [1].
#
# [1] Colas Le Guernic: Reachability analysis of hybrid systems with linear
#     continuous dynamics. PhD thesis, 2009.
# ==============================================================================
using MathematicalSystems, LazySets

function fiveDimSys_model()
    # system matrix
    D = [-1.0 -4.0 0.0 0.0 0.0;
         4.0 -1.0 0.0 0.0 0.0;
         0.0 0.0 -3.0 1.0 0.0;
         0.0 0.0 -1.0 -3.0 0.0;
         0.0 0.0 0.0 0.0 -2.0]
    P = [0.6 -0.1 0.1 0.7 -0.2;
         -0.5 0.7 -0.1 -0.8 0.0;
         0.9 -0.5 0.3 -0.6 0.1;
         0.5 -0.7 0.5 0.6 0.3;
         0.8 0.7 0.6 -0.3 0.2]
    A = P * D * inv(P)

    # state domain
    X = Universe(5)

    # input domain
    U = Ball2(zeros(5), 0.01)

    # continuous LTI system
    S = @system(x' = Ax + u, x ∈ X, u ∈ U)

    return S
end
