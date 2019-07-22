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
    dx[2] = a*b*x[1] - a*x[2]
end

function spiking_neurons()

    automaton = LightAutomaton(1) # one node

    inv_one = HPolyhedron([HalfSpace([1.0, 0.0], 30.0)])


    m1 = ConstrainedBlackBoxContinuousSystem(spiking_neuron!, 2, inv_one)

    modes = [m1]

    add_transition!(automaton, 1, 1, 1)

    guard_alpha = HPolyhedron([HalfSpace([-1.0, 0.0], -30.0)])


    A = [0.0 0.0;0.0 1.0]
    b = [-65.0, 8.0]
    t1 = ConstrainedAffineMap(A, b, guard_alpha)
    #resetmaps
    resetmaps = [t1]
    # switching
    switching = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switching)
    #initial condition in mode one
    X0 = Hyperrectangle(low = [-65.0, -0.2], high = [-60.0, 0.2])

    initial_condition = [(1, X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    options = Options(:mode=>"reach", :T=>100.0, :plot_vars=>[1, 2], :project_reachset=>false)

    return (system, options)
end

problem, options = spiking_neurons();

options = Options(:mode=>"reach", :T=>100.0, :plot_vars=>[1, 2], :project_reachset=>false, :verbosity => "info")

opC = TMJets(:orderT=>5, :orderQ=>2, :abs_tol=>1e-10)
opD = LazyDiscretePost(:check_invariant_intersection=>true)
@time sol_TMJets = solve(problem, options, opC, opD)

using IntervalArithmetic
a = (100/0.84)*IntervalArithmetic.Interval(sol_TMJets.Xk[1].t_start,sol_TMJets.Xk[1].t_end)
b = IntervalArithmetic.Interval(sol_TMJets.Xk[1].X.radius[1],sol_TMJets.Xk[1].X.radius[2])
fig1 = plot(a×b, colour = "green")
for i =2:length(sol_TMJets.Xk)
    a = (100/0.84)*IntervalArithmetic.Interval(sol_TMJets.Xk[i].t_start,sol_TMJets.Xk[i].t_end)
    b = IntervalArithmetic.Interval(sol_TMJets.Xk[i].X.radius[1],sol_TMJets.Xk[i].X.radius[2])
    fig1 =  plot!(a×b, colour = "green")
end

a = (100/0.84)*IntervalArithmetic.Interval(sol_TMJets.Xk[1].t_start,sol_TMJets.Xk[1].t_end)
b = IntervalArithmetic.Interval(sol_TMJets.Xk[1].X.center[1],sol_TMJets.Xk[1].X.center[2])
fig2 = plot(a×b, colour = "green")
for i =2:length(sol_TMJets.Xk)
    a = (100/0.84)*IntervalArithmetic.Interval(sol_TMJets.Xk[i].t_start,sol_TMJets.Xk[i].t_end)
    b = IntervalArithmetic.Interval(sol_TMJets.Xk[i].X.center[1],sol_TMJets.Xk[i].X.center[2])
    fig2 =  plot!(a×b, colour = "green")
end
