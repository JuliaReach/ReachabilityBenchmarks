# =============================================================================
# Spiking neuron
#
# system type: hybrid system with polynomial dynamics
# state dimension: 2
# modes: 1
# transitions: 1
#
# This model is adapted from Izhikevich's model [1].
# In [2, 3] this model has been formalized as a hybrid automaton.
# 
# [1] Eugene M. Izhikevich: Dynamical systems in neuroscience. MIT press, 2007.
# [2] https://flowstar.org/examples/
# [3] Xin Chen: Reachability analysis of non-linear hybrid systems using Taylor
# models. PhD thesis, 2015.
# =============================================================================
using HybridSystems, MathematicalSystems, LazySets, TaylorIntegration

@taylorize function flow!(dx, x, params, t)
    local a = 0.02
    local b = 0.2
    local I = 40
    dx[1] = (0.04*(x[1]*x[1]) + 5*x[1]) + ((140+I)-x[2])
    dx[2] = a*((b*x[1]) - x[2])
    return dx
end

function spikingNeuron_model()
    # hybrid automaton
    HA = LightAutomaton(1)

    # mode 1
    X = HPolyhedron([HalfSpace([1.0, 0.0], 30.0)])  # x1 ≤ 30
    m1 = @system(x' = flow!(x), dim: 2, x ∈ X)

    # transition mode 1 → mode 1 (self loop)
    add_transition!(HA, 1, 1, 1)
    G = HPolyhedron([HalfSpace([-1.0, 0.0], -30.0)])  # x1 ≥ 30
    A = [0.0 0.0; 0.0 1.0]
    b = [-65.0, 8.0]
    R11 = ConstrainedAffineMap(A, b, G)  # x1 := -65, x2 := x2 + 8

    # hybrid system
    S = HybridSystem(HA, [m1], [R11], [AutonomousSwitching()])

    return S
end
