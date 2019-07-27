# ==============================================================================

# ==============================================================================

using Reachability, HybridSystems
using Plots

@taylorize function transmission_line_circuits_one!(t, x, dx)
   local α = 5
   dx[1] = -2*x[1] + x[2] + 2 - exp(α*x[1]) - exp(α*(x[1] - x[2])) + x[7]
   dx[2] = -2*x[2] + x[1] + x[3] + exp(α*(x[1] - x[2])) - exp(α*(x[2] - x[3])) 
   dx[3] = -2*x[3] + x[2] + x[4] + exp(α*(x[2] - x[3])) - exp(α*(x[3] - x[4]))
   dx[4] = -2*x[4] + x[3] + x[5] + exp(α*(x[3] - x[4])) - exp(α*(x[4] - x[5]))
   dx[5] = -2*x[5] + x[4] + x[6] + exp(α*(x[4] - x[5])) - exp(α*(x[5] - x[6]))
   dx[6] = -x[6] + x[5] -1 + exp(α*(x[5] - x[6]))
   return dx
end

@taylorize function transmission_line_circuits_two!(t, x, dx)
   local α = 5
   dx[1] = -2*x[1] + x[2] + 2 - exp(α*x[1]) - exp(α*(x[1] - x[2])) + x[7]
   dx[2] = -2*x[2] + x[1] + x[3] + exp(α*(x[1] - x[2])) - exp(α*(x[2] - x[3])) 
   dx[3] = -2*x[3] + x[2] + x[4] + exp(α*(x[2] - x[3])) - exp(α*(x[3] - x[4]))
   dx[4] = -2*x[4] + x[3] + x[5] + exp(α*(x[3] - x[4])) - exp(α*(x[4] - x[5]))
   dx[5] = -2*x[5] + x[4] + x[6] + exp(α*(x[4] - x[5])) - exp(α*(x[5] - x[6]))
   dx[6] = -x[6] + x[5] -1 + exp(α*(x[5] - x[6]))
   return dx
end

@taylorize function transmission_line_circuits_three!(t, x, dx)
   local α = 5
   dx[1] = -2*x[1] + x[2] + 2 - exp(α*x[1]) - exp(α*(x[1] - x[2])) + x[7]
   dx[2] = -2*x[2] + x[1] + x[3] + exp(α*(x[1] - x[2])) - exp(α*(x[2] - x[3])) 
   dx[3] = -2*x[3] + x[2] + x[4] + exp(α*(x[2] - x[3])) - exp(α*(x[3] - x[4]))
   dx[4] = -2*x[4] + x[3] + x[5] + exp(α*(x[3] - x[4])) - exp(α*(x[4] - x[5]))
   dx[5] = -2*x[5] + x[4] + x[6] + exp(α*(x[4] - x[5])) - exp(α*(x[5] - x[6]))
   dx[6] = -x[6] + x[5] -1 + exp(α*(x[5] - x[6]))
   return dx
end

function circuits()

   automaton = LightAutomaton(2) # two modes

   inv_one = HPolyhedron([HalfSpace([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 1.0)]) 

   inv_two = HPolyhedron([HalfSpace([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 2.0), 
                          HalfSpace([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0], -1.0)]) 

   inv_three = HPolyhedron([HalfSpace([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0], -2.0)]) 
                                                                                                                                          

   m1 = ConstrainedBlackBoxContinuousSystem(transmission_line_circuits_one!, 3, inv_one)

   m2 = ConstrainedBlackBoxContinuousSystem(transmission_line_circuits_two!, 3, inv_two)

   m3 = ConstrainedBlackBoxContinuousSystem(transmission_line_circuits_three!, 3, inv_three)

   modes = [m1, m2, m3]  #modes


   add_transition!(automaton, 1, 2, 1)    # first transion
   add_transition!(automaton, 2, 3, 2)    # second transion
   add_transition!(automaton, 2, 1, 3)    # third transion
   add_transition!(automaton, 3, 2, 4)    # fourth transion 
     
   G_first = HPolyhedron([HalfSpace([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0], -1.0)]) 
   G_second = HPolyhedron([HalfSpace([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0], -2.0)]) 
   G_third = HPolyhedron([HalfSpace([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 1.0)]) 
   G_fourth = HPolyhedron([HalfSpace([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 2.0)]) 

   A1 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0]
   A2 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
        
   
   b1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0]
   b2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0]
   b3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
   
   R_first = ConstrainedAffineMap(A1, b2, G_first)
   R_second = ConstrainedAffineMap(A2, b3, G_second)
   R_third = ConstrainedAffineMap(A2, b1, G_third)
   R_fourth = ConstrainedAffineMap(A1, b2, G_fourth)
   # resetmaps
   resetmaps = [R_first, R_second, R_third, R_fourth]
   # switchings
   switching = AutonomousSwitching()
   switchings = fill(switching, 4)

   ℋ = HybridSystem(automaton, modes, resetmaps, switchings)
   # initial condition in mode three
   X0 = Hyperrectangle(low = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], high = [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01])
   # initial condition 
   initial_condition = [(1, X0)]  

   system = InitialValueProblem(ℋ, initial_condition)

   options = Options(:mode=>"reach", :T=>3.5, :plot_vars=>[1, 2], :project_reachset=>false)

   return (system, options)
end


problem, options = circuits();

options = Options(:mode=>"reach", :T=>3.0, :plot_vars=>[1, 2], :project_reachset=>false, :verbosity => "info")

@time sol_TMJets = solve(problem, options, TMJets(:orderT=>5, :orderQ=>2, :abs_tol=>1e-10),LazyDiscretePost(:check_invariant_intersection=>true));

plot(sol_TMJets, xlab="x", ylab="v", alpha=.5, color=:lightblue)
