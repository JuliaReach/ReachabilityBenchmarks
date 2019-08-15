using Reachability, HybridSystems, MathematicalSystems, LazySets, LinearAlgebra, SX, SymEngine
using Plots, LaTeXStrings

file = "models/LinearSwitcher/model.xml"
model = readsxmodel(file, raw_dict=true)

function linear_switching(; X0 = Singleton([3.1, 4.0, 0.0, 0.0, 0.0]),
                            U = Interval(-1.0, 1.0),
                            T = 1.0,
                            ε = 1e-6)

    n = 5 # number of state variables
    m = 1 # number of input variables
    state_vars = convert.(Basic, [fi.args[1].args[1] for fi in model["flows"][1]])
    input_vars = [convert(Basic, collect(keys(model["variables"]))[6])] # same as convert(Basic, :u)

    id_location = 1
    A, B, c = get_coeffs(model["flows"][id_location], n, m, state_vars, input_vars) # new function in SX
    X = HPolyhedron([HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -3.0 + ε)]) # x >= 3
    q1 = ConstrainedLinearControlContinuousSystem(A, B, X, U)

    id_location = 2
    A, B, c = get_coeffs(model["flows"][id_location], n, m, state_vars, input_vars)
    X = HPolyhedron([HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -2.0 + ε)]) # x >= 2
    q2 = ConstrainedLinearControlContinuousSystem(A, B, X, U)

    id_location = 3
    A, B, c = get_coeffs(model["flows"][id_location], n, m, state_vars, input_vars)
    X = HPolyhedron([HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -1.0 + ε)]) # x >= 1
    q3 = ConstrainedLinearControlContinuousSystem(A, B, X, U)

    id_location = 4
    A, B, c = get_coeffs(model["flows"][id_location], n, m, state_vars, input_vars) # new function in SX
    X = HPolyhedron([HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], 0.0 + ε)]) # x >= 0
    q4 = ConstrainedLinearControlContinuousSystem(A, B, X, U)

    id_location = 5
    A, B, c = get_coeffs(model["flows"][id_location], n, m, state_vars, input_vars) # new function in SX
    X = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 1.0+ε)]) # x <= 1
    q5 = ConstrainedLinearControlContinuousSystem(A, B, X, U)

    # automaton structure
    automaton = LightAutomaton(5)

    modes = [q1, q2, q3, q4, q5]

    # transitions
    add_transition!(automaton, 1, 2, 1)
    add_transition!(automaton, 2, 3, 2)
    add_transition!(automaton, 3, 4, 3)
    add_transition!(automaton, 4, 5, 4)
    add_transition!(automaton, 5, 1, 5)

    # guards
    G12 = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 3.0 + ε),
                       HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -3.0 + ε)]) # x1 = 3

    G23 = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 2.0 + ε),
                       HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -2.0 + ε)]) # x1 = 2

    G34 = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 1.0 + ε),
                       HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -1.0 + ε)]) # x1 = 1

    G45 = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 0.0 + ε),
                       HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -0.0 + ε)]) # x1 = 0

    G51 = HPolyhedron([HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], 1.0 + ε),
                       HalfSpace([-1.0, 0.0, 0.0, 0.0, 0.0], -1.0 + ε)]) # x1 = 1

    resetmaps = [ConstrainedIdentityMap(2, G12),
                 ConstrainedIdentityMap(2, G23),
                 ConstrainedIdentityMap(2, G34),
                 ConstrainedIdentityMap(2, G45),
                 ConstrainedIdentityMap(2, G51)];

    # switching
    switchings = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition in mode q1
    initial_condition = [(1, X0)]

    problem = InitialValueProblem(ℋ, initial_condition)

    border = HPolyhedron([HalfSpace([1.0; zeros(4)], -1.2)]) # x1 <= -1.2

    property = BadStatesProperty(border)

    options = Options(:mode=>"reach", :T=>T, :ε_proj=>1e-3,:set_type_proj=>HPolygon, :overapproximation=>Approximations.OctDirections, :plot_vars=>[1, 2], :project_reachset=>false, :property=>property)

    return problem, options
end

X0 = Singleton([3.1, 4.0, 0.0, 0.0, 0.0])
U = Interval(-1.0, 1.0)
problem, options = linear_switching(; X0=X0, U=U, T=1.0, ε=1e-6);

using LazySets.Approximations
using LazySets.Approximations: project, overapproximate

@time begin
    opC = BFFPS19(:δ=>0.001, :partition=>[1:2, 3:3, 4:4,5:5],:ε_proj=>0.001)
    opD = DecomposedDiscretePost(:out_vars=>[1,2], :clustering=>:none)
    sol = solve(problem, options, opC, opD)
end;
