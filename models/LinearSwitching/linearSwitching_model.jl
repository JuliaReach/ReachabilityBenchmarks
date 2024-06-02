# ==============================================================================
# Random linear switching system
#
# system type: hybrid system with linear dynamics
# state dimension: 5
# input dimension: 1
# modes: 5
# transitions: 5
#
# The model is taken from [1].
#
# [1] https://ths.rwth-aachen.de/research/projects/hypro/5-dimensional-switching-linear-system/
# ==============================================================================
using HybridSystems, MathematicalSystems, LazySets, SymEngine
using SX: readsxmodel, _get_coeffs

function linearSwitching_model()
    n = 5  # state dimension
    U = Interval(-1.0, 1.0)  # common input domain
    m = dim(U)  # input dimension
    ε = 1e-6  # auxiliary bloating of guards for ensuring intersection

    # hybrid automaton
    HA = GraphAutomaton(5)

    # SpaceEx model
    file = @current_path "SpaceEx/model.xml"
    model = readsxmodel(file; raw_dict=true)
    variables = convert.(Basic, [f.args[1].args[1] for f in model["flows"][1]])
    inputs = [convert(Basic, :u)]
    load_dynamics(loc) = _get_coeffs(model["flows"][loc], n, m, variables, inputs)

    # mode 1
    A, B, _ = load_dynamics(1)
    X = HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -3.0 + ε)  # x1 ≥ 3
    q1 = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    # mode 2
    A, B, c = load_dynamics(2)
    X = HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -2.0 + ε)  # x1 ≥ 2
    q2 = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    # mode 3
    A, B, _ = load_dynamics(3)
    X = HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -1.0 + ε)  # x1 ≥ 1
    q3 = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    # mode 4
    A, B, _ = load_dynamics(4)
    X = HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], 0.0 + ε)  # x1 ≥ 0
    q4 = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    # mode 5
    A, B, _ = load_dynamics(5)
    X = HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 1.0 + ε)  # x1 ≤ 1
    q5 = @system(x' = Ax + Bu, x ∈ X, u ∈ U)

    # transition mode 1 → mode 2
    add_transition!(HA, 1, 2, 1)
    G = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 3.0 + ε),
                     HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -3.0 + ε)])  # x1 = 3
    R12 = ConstrainedIdentityMap(5, G)

    # transition mode 2 → mode 3
    add_transition!(HA, 2, 3, 2)
    G = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 2.0 + ε),
                     HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -2.0 + ε)])  # x1 = 2
    R23 = ConstrainedIdentityMap(5, G)

    # transition mode 3 → mode 4
    add_transition!(HA, 3, 4, 3)
    G = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 1.0 + ε),
                     HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -1.0 + ε)])  # x1 = 1
    R34 = ConstrainedIdentityMap(5, G)

    # transition mode 4 → mode 5
    add_transition!(HA, 4, 5, 4)
    G = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 0.0 + ε),
                     HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -0.0 + ε)])  # x1 = 0
    R45 = ConstrainedIdentityMap(5, G)

    # transition mode 5 → mode 1
    add_transition!(HA, 5, 1, 5)
    G = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 1.0 + ε),
                     HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -1.0 + ε)])  # x1 = 1
    R51 = ConstrainedIdentityMap(5, G)

    # hybrid system
    S = HybridSystem(HA, [q1, q2, q3, q4, q5], [R12, R23, R34, R45, R51],
                     fill(AutonomousSwitching(), 5))

    return S
end
