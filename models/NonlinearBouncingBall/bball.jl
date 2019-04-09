# ==============================================================================
# Bouncing ball model with air friction
#
# See Example 4.1.2 page 98 in Chen, X. (2015).
# Reachability analysis of non-linear hybrid systems using taylor models
# (Doctoral dissertation, Fachgruppe Informatik, RWTH Aachen University).
# https://www.cs.colorado.edu/~xich8622/papers/thesis.pdf
# ==============================================================================

using Revise, MathematicalSystems, Reachability, LinearAlgebra, HybridSystems
using LazySets
using Reachability: solve
using Plots
using TaylorIntegration

@taylorize function bball_up!(t, x, dx)
    dx[1] = x[2]
    dx[2] = -9.8 - 0.1*(x[1])^2
    return dx
end

@taylorize function bball_down!(t, x, dx)
    dx[1] = x[2]
    dx[2] = -9.8 + 0.1*(x[1])^2
    return dx
end

function bouncing_ball()

    automaton = LightAutomaton(2) # two modes

    inv_up = HPolyhedron([HalfSpace([-1.0, 0.0], 0.0),  # x >= 0
                          HalfSpace([0.0, -1.0], 0.0)]) # v >= 0

    inv_down = HPolyhedron([HalfSpace([-1.0, 0.0], 0.0),  # x >= 0
                            HalfSpace([ 0.0, 1.0], 0.0)])  # v <= 0

    m1 = ConstrainedBlackBoxContinuousSystem(bball_up!, 2, inv_up)

    m2 = ConstrainedBlackBoxContinuousSystem(bball_down!, 2, inv_down)

    modes = [m1, m2]  #modes

    add_transition!(automaton, 2, 1, 1)    # alpha transition
    add_transition!(automaton, 1, 2, 2)    # beta transition

    Gα = HPolyhedron([HalfSpace([1.0, 0.0], 0.0),    # x>=0
                      HalfSpace([-1.0, 0.0], 0.0)])  # x<=0
    Gβ =  HPolyhedron([HalfSpace([0.0, 1.0], 0.0),   # v<=0
                       HalfSpace([0.0,-1.0], 0.0)])  # v>=0


    A = [1.0 0.0; 0.0 -0.8]
    Rα = ConstrainedLinearMap(A, Gα)
    Rβ = ConstrainedIdentityMap(2, Gβ)

    # resetmaps
    resetmaps = [Rα, Rβ]

    # switchings
    switching = AutonomousSwitching()
    switchings = fill(switching, 2)

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition in mode two
    X0 = Hyperrectangle(low=[4.9, -0.2], high=[5.1, 0.0])

    initial_condition = [(2, X0)]  # initial condition in the "down" mode

    system = InitialValueProblem(ℋ, initial_condition)

    options = Options(:mode=>"reach", :T=>6.0, :plot_vars=>[1, 2], :project_reachset=>false)

    return (system, options)
end


problem, options = bouncing_ball();

options = Options(:mode=>"reach", :T=>3.0, :plot_vars=>[1, 2], :project_reachset=>false, :verbosity => "info")

@time sol_TMJets = solve(problem, options, TMJets(:orderT=>5, :orderQ=>2, :abs_tol=>1e-10),LazyDiscretePost(:check_invariant_intersection=>true));

plot(sol_TMJets, xlab="x", ylab="v", alpha=.5, color=:lightblue)
