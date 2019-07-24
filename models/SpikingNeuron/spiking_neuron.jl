#=
This is a spiking neuron model adapted from Izhikevich's model [1].
In [2, 3] this model has been formalized as a hybrid automaton consisting of
one mode and two variables.

[1] Izhikevich, Eugene M. Dynamical systems in neuroscience. MIT press, 2007.
[2] https://flowstar.org/examples/
[3] Chen, Xin. Reachability Analysis of Non-Linear Hybrid Systems Using Taylor Models.
    Diss. Fachgruppe Informatik, RWTH Aachen University, 2015. 
=#
using MathematicalSystems, Reachability, HybridSystems, Plots

@taylorize function spiking_neuron!(t, x, dx)
    local a = 0.02
    local b = 0.2
    local I = 40
    dx[1] = (0.04*(x[1]*x[1]) + 5*x[1]) + ((140+I)-x[2])
    dx[2] = a*((b*x[1]) - x[2])
end

function spiking_neuron()
    # mode and invariant
    automaton = LightAutomaton(1)
    invariant = HPolyhedron([HalfSpace([1.0, 0.0], 30.0)])
    modes = [ConstrainedBlackBoxContinuousSystem(spiking_neuron!, 2, invariant)]

    # transition and guard
    add_transition!(automaton, 1, 1, 1)
    guard = HPolyhedron([HalfSpace([-1.0, 0.0], -30.0)])

    # reset map and switching
    A = [0.0 0.0; 0.0 1.0]
    b = [-65.0, 8.0]
    resetmaps = [ConstrainedAffineMap(A, b, guard_alpha)]
    switching = [AutonomousSwitching()]

    # create hybrid system
    ℋ = HybridSystem(automaton, modes, resetmaps, switching)

    #initial condition
    X0 = Hyperrectangle(low=[-65.0, -0.2], high=[-60.0, 0.2])
    initial_condition = [(1, X0)]

    system = InitialValueProblem(ℋ, initial_condition)
    options = Options(:mode=>"reach", :T=>100.0, :plot_vars=>[1, 2], :project_reachset=>false)

    return (system, options)
end

problem, options = spiking_neuron()

opC = TMJets(:orderT=>5, :orderQ=>2, :abs_tol=>1e-10)
opD = LazyDiscretePost(:check_invariant_intersection=>true)
@time sol_TMJets = solve(problem, options, opC, opD)

#=
using IntervalArithmetic
a = (100/0.84)*IntervalArithmetic.Interval(sol_TMJets.Xk[1].t_start, sol_TMJets.Xk[1].t_end)
b = IntervalArithmetic.Interval(sol_TMJets.Xk[1].X.radius[1], sol_TMJets.Xk[1].X.radius[2])
fig1 = plot(a×b, colour = "green")
for i =2:length(sol_TMJets.Xk)
    a = (100/0.84)*IntervalArithmetic.Interval(sol_TMJets.Xk[i].t_start, sol_TMJets.Xk[i].t_end)
    b = IntervalArithmetic.Interval(sol_TMJets.Xk[i].X.radius[1], sol_TMJets.Xk[i].X.radius[2])
    fig1 =  plot!(a×b, colour = "green")
end

a = (100/0.84)*IntervalArithmetic.Interval(sol_TMJets.Xk[1].t_start, sol_TMJets.Xk[1].t_end)
b = IntervalArithmetic.Interval(sol_TMJets.Xk[1].X.center[1], sol_TMJets.Xk[1].X.center[2])
fig2 = plot(a×b, colour = "green")
for i =2:length(sol_TMJets.Xk)
    a = (100/0.84)*IntervalArithmetic.Interval(sol_TMJets.Xk[i].t_start, sol_TMJets.Xk[i].t_end)
    b = IntervalArithmetic.Interval(sol_TMJets.Xk[i].X.center[1], sol_TMJets.Xk[i].X.center[2])
    fig2 =  plot!(a×b, colour = "green")
end
=#
