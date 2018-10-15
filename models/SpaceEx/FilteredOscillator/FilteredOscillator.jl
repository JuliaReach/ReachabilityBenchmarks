# Filtered Oscillator
# See: https://flowstar.org/benchmarks/filtered-oscillator/
# ============================

using HybridSystems, MathematicalSystems, LazySets, Reachability, Polyhedra
import LazySets.HalfSpace


n0 = 2
system_dimension = n0 + 2
z = zeros(n0)

# Transition graph (automaton)
a = LightAutomaton(4);
add_transition!(a, 3, 4, 1);
add_transition!(a, 4, 2, 2);
add_transition!(a, 2, 1, 3);
add_transition!(a, 1, 3, 4);

# common flow
A = zeros(system_dimension, system_dimension)
A[1,1], A[2,2] = -2., -1.
A[3,1], A[3,3] = 5., -5.
for i = 4 : system_dimension
    A[i,i-1], A[i,i] = 5., -5.
end

# common U
U = Singleton([1.0]);
# Modes

#Mode 1
B = [1.4; -0.7; z]
X = HPolyhedron([HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
               HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
m_1 = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

#Mode 2
B = [-1.4; 0.7; z]
X = HPolyhedron([HalfSpace([1.0; 0.0; z], 0.0),  # x <= 0
               HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
m_2 = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

#Mode 3
B = [1.4; -0.7; z]
X = HPolyhedron([HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
               HalfSpace([-0.714286; -1.0; z], 0.0)])  # 0.714286*x + y >= 0
m_3 = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

#Mode 4
B = [-1.4; 0.7; z]
X = HPolyhedron([HalfSpace([0.714286; 1.0; z], 0.0),  # 0.714286*x + y <= 0
               HalfSpace([-1.0; 0.0; z], 0.0)])  # x >= 0
m_4 = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

m = [m_1, m_2, m_3, m_4];

#Transitions

# common resets
A_trans = eye(system_dimension)

# Transition l3 -> l4
X_l3l4 = HPolyhedron([HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
               HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
               HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
r1 = ConstrainedLinearDiscreteSystem(A_trans, X_l3l4);
# Transition l4 -> l2
X_l4l2 = HPolyhedron([HalfSpace([0.714286; 1.0; z], 0.0),  # 0.714286*x + y <= 0
               HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
               HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
r2 = ConstrainedLinearDiscreteSystem(A_trans, X_l4l2);
# Transition l2 -> l1
X_l2l1 = HPolyhedron([HalfSpace([1.0; 0.0; z], 0.0),  # x <= 0
               HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
               HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
r3 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l1);
# Transition l1 -> l3
X_l1l3 = HPolyhedron([HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
               HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
               HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
r4 = ConstrainedLinearDiscreteSystem(A_trans, X_l1l3);

r = [r1,r2,r3,r4]
# Switchings
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s);

# initial condition in mode 1
X0 = Hyperrectangle(low=[0.2; -0.1 * ones(system_dimension-1)],
                    high=[0.4; 0.1 * ones(system_dimension-1)]);

system = InitialValueProblem(HS, [(3, X0)]);
plot_vars = [1, 2]
options_common = Options(:mode=>"reach",:vars=>1:system_dimension, :T=>10.0, :δ=>0.01,
                  :plot_vars=>plot_vars, :max_jumps=>4,
                  :partition=>[1:system_dimension], :ε_proj=>0.001,
                  :clustering=>:chull, :verbosity=>1);

options_proj_false = Options(:project_reachset=>false)
options_proj_true = Options(:project_reachset=>true)

options_pr_t = merge!(options_common, options_proj_true)
options_pr_f = merge!(options_common, options_proj_false)

# default algorithm
sol = solve(system, options_pr_t);
# specify lazy discrete post operator
sol = solve(system, options_pr_f, Reachability.BFFPSV18(),
            Reachability.ReachSets.LazyTextbookDiscretePost());

N = Float64
sol_processed =  Reachability.ReachSolution([Reachability.ReachSet{CartesianProductArray{N}, N}(
            CartesianProductArray{N, HPolytope{N}}([LazySets.Approximations.overapproximate(rs.X, LazySets.Approximations.OctDirections(4))]),
            rs.t_start, rs.t_end) for rs in sol.Xk], sol.options)

sol_proj = Reachability.ReachSolution(Reachability.project_reach(
    sol_processed.Xk, plot_vars, system_dimension, sol.options), sol.options);

# specify overapproximating discrete post operator
sol = solve(system, options_pr_t, Reachability.BFFPSV18(),
           Reachability.ReachSets.ApproximatingDiscretePost(
                Options(:overapproximation=>Hyperrectangle)));
