# =============================================================================
# Bouncing ball with air friction
#
# system type: hybrid system with polynomial dynamics
# state dimension: 2
# modes: 2
# transitions: 2
#
# See Example 4.1.2 (page 98f) in [1].
# 
# [1] Xin Chen: Reachability analysis of non-linear hybrid systems using Taylor
# models. PhD thesis, 2015.
# =============================================================================
using HybridSystems, MathematicalSystems, LazySets, TaylorIntegration

@taylorize function flow_down!(dx, x, params, t)
    dx[1] = x[2]
    dx[2] = -9.8 + 0.1*(x[2])^2
    return dx
end

@taylorize function flow_up!(dx, x, params, t)
    dx[1] = x[2]
    dx[2] = -9.8 - 0.1*(x[2])^2
    return dx
end

function bouncingBallNonlinear_model()
    # hybrid automaton with state variables x, v
    HA = LightAutomaton(2)

    # mode 1 ("down")
    X = HPolyhedron([HalfSpace([-1.0, 0.0], 0.0),   # x ≥ 0
                     HalfSpace([ 0.0, 1.0], 0.0)])  # v ≤ 0
    m1 = @system(x' = flow_down!(x), dim: 2, x ∈ X)

    # mode 2 ("up")
    X = HPolyhedron([HalfSpace([-1.0, 0.0], 0.0),   # x ≥ 0
                     HalfSpace([0.0, -1.0], 0.0)])  # v ≥ 0
    m2 = @system(x' = flow_up!(x), dim: 2, x ∈ X)

    # α transition down → up
    add_transition!(HA, 1, 2, 1)
    G = Hyperplane([1.0, 0.0], 0.0)  # x = 0
    A = [1.0 0.0; 0.0 -0.8]
    Rα = ConstrainedLinearMap(A, G)  # v := -0.8v

    # β transition up → down
    add_transition!(HA, 2, 1, 2)
    G = Hyperplane([0.0, 1.0], 0.0)  # v = 0
    Rβ = ConstrainedIdentityMap(2, G)

    # hybrid system
    S = HybridSystem(HA, [m1, m2], [Rα, Rβ], fill(AutonomousSwitching(), 2))

    return S
end
