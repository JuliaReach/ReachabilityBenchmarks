# Filtered Oscillator
# See: https://flowstar.org/benchmarks/filtered-oscillator/
# ============================

using HybridSystems, MathematicalSystems, LazySets, Polyhedra
import LazySets.HalfSpace


system_dimension = 32 + 2
# Transition graph (automaton)
a = LightAutomaton(4);
add_transition!(a, 3, 4, 1);
add_transition!(a, 4, 2, 1);
add_transition!(a, 2, 1, 1);
add_transition!(a, 1, 3, 1);


flow = zeros(eye(Float64, system_dimension))

flow[3,1], flow[3,3] = 5.,-5.
for i = 4 : system_dimension
    flow[i,i-1], flow[i,i] = 5.,-5.
end

# Modes

#Mode 1
flow[1,1], flow[2,2] = -2., -1.
A = flow
B = reshape(vcat([1.4, -0.7], zeros(32)), (34, 1))
X = HalfSpace([-1.0, 0.0], 0.0); # x >= 0
U = Singleton([1.0]);
m_1 = [ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U)];

#Mode 2
flow[1,1], flow[2,2] = -2., -1.
A = flow
B = reshape(vcat([-1.4, 0.7], zeros(32)), (34, 1))
X = HalfSpace([-1.0, 0.0], 0.0); # x >= 0
U = Singleton([1.0]);
m_2 = [ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U)];

#Mode 3
flow[1,1], flow[2,2] = -2., -1.
A = flow
B = reshape(vcat([1.4, -0.7], zeros(32)), (34, 1))
X = HalfSpace([-1.0, 0.0], 0.0); # x >= 0
U = Singleton([1.0]);
m_3 = [ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U)];

#Mode 4
flow[1,1], flow[2,2] = -2., -1.
A = flow
B = reshape(vcat([-1.4, 0.7], zeros(32)), (34, 1))
X = HalfSpace([-1.0, 0.0], 0.0); # x >= 0
U = Singleton([1.0]);
m_4 = [ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U)];

m = [m_1, m_2, m_3, m_4];

#Transitions

# Transition l3 -> l4
A_l3l4 = ones(eye(Float64, system_dimension))
X_l3l4 = HPolytope([HalfSpace(vcat([-1.0, 0.0], zeros(32)), 0.0),   # x >= 0
               HalfSpace(vcat([-0.714286, -1.0], zeros(32)), 0.0),  # 0.714286*x + y  >= 0
               HalfSpace(vcat([0.714286, 0.0], zeros(32)), 0.0)])  # 0.714286*x + y <= 0
r1 = [ConstrainedLinearDiscreteSystem(A_l3l4, X_l3l4)];
# Transition l4 -> l2
A_l4l2 = ones(eye(Float64, system_dimension))
X_l4l2 = HPolytope([HalfSpace(vcat([0.714286, 0.0], zeros(32)), 0.0),   # 0.714286*x + y <= 0
               HalfSpace(vcat([-1.0, 0.0], zeros(32)), 0.0),  # x >= 0
               HalfSpace(vcat([1.0, 0.0], zeros(32)), 0.0)]) # x <= 0
r2 = [ConstrainedLinearDiscreteSystem(A_l4l2, X_l4l2)];
# Transition l2 -> l1
A_l2l1 = ones(eye(Float64, system_dimension))
X_l2l1 = HPolytope([HalfSpace(vcat([1.0, 0.0], zeros(32)), 0.0),   # x <= 0
               HalfSpace(vcat([-0.714286, -1.0], zeros(32)), 0.0),  #  0.714286*x + y >= 0
               HalfSpace(vcat([0.714286, 0.0], zeros(32)), 0.0)]) #  0.714286*x + y <= 0
r3 = [ConstrainedLinearDiscreteSystem(A_l2l1, X_l2l1)];
# Transition l1 -> l3
A_l1l3 = ones(eye(Float64, system_dimension))
X_l1l3 = HPolytope([HalfSpace(vcat([-0.714286, -1.0], zeros(32)), 0.0),   # 0.714286*x + y >= 0
               HalfSpace(vcat([-1.0, 0.0], zeros(32)), 0.0),  # x >= 0
               HalfSpace(vcat([1.0, 0.0], zeros(32)), 0.0)]) # x <= 0
r4 = [ConstrainedLinearDiscreteSystem(A_l1l3, X_l1l3)];

r = [r1,r2,r3,r4]
# Switchings
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s);

# initial condition in mode 1
X0 = Hyperrectangle(low=vcat([0.2, -0.1], zeros(32)), high=vcat([0.3, 0.1], zeros(32)));

system = InitialValueProblem(HS, X0);
input_options = Options(:mode=>"reach");

problem_options = Options(:vars=>[1,2], :T=>15.0, :δ=>0.05, :plot_vars=>[1, 2],
                          :max_jumps=>20, :verbosity=>1);
options_input = merge(problem_options, input_options);
sol = solve(system, options_input);
