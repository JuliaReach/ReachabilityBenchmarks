using Revise, MathematicalSystems,Reachability, LinearAlgebra, HybridSystems
using LazySets;
using Reachability: solve
using Plots
using TaylorIntegration

@taylorize function spiking_neuron!(t, x, dx)
    local a = 0.02
    local b = 0.2
    local I = 40
    dx[1] = 0.04*x[1]^2 + 5*x[1] + 140 - x[2] + I
    dx[2] = a*(b*x[1] - x[2])
end

function spiking_neurons()

    automaton = LightAutomaton(1) # one node

    inv_one = HPolyhedron([HalfSpace([-1.0, 0.0], -0.2),    #x[1] >= -0.2
                           HalfSpace([0.0, -1.0], -70.0)])  #x[2] >= -70

    m1 = ConstrainedBlackBoxContinuousSystem(spiking_neuron!, 2, inv_one)

    modes = [m1]

    add_transition!(automaton, 1, 1, 1)

    guard_alpha = HPolyhedron([HalfSpace([-1.0, 0.0], 30.0),   # x[1] >= 30
                               HalfSpace([1.0, 0.0], 30.0)]) # x[1] <= 30

    A = [1.0 0.0;0.0 1.0]
    b = [-65.0, 8.0]
    t1 = ConstrainedAffineMap(A,b,guard_alpha)
    #resetmaps
    resetmaps = [t1]
    # switching
    switching = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switching)
    #initial condition in mode one
    X0 = Hyperrectangle(low = [-0.2, -65.0], high = [0.2, -60.0])

    initial_condition = [(1, X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    options = Options(:mode=>"reach", :T=>100.0, :plot_vars=>[1, 2], :project_reachset=>false)

    return (system, options)
end

problem, options = spiking_neurons();

options = Options(:mode=>"reach", :T=>100.0, :plot_vars=>[1, 2], :project_reachset=>false, :verbosity => "info")

@time sol_TMJets = solve(problem, options, TMJets(:orderT=>5, :orderQ=>2, :abs_tol=>1e-10),LazyDiscretePost(:check_invariant_intersection=>true))

plot(sol_TMJets, use_subindices=false, alpha=.5, color=:lightblue)
