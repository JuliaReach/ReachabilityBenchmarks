# ==============================================================================
# Filtered Oscillator
#
# See: https://flowstar.org/benchmarks/filtered-oscillator/
# ==============================================================================

using HybridSystems, MathematicalSystems, LazySets, Reachability, Polyhedra, Optim
using LinearAlgebra, SparseArrays

import LazySets.HalfSpace
import LazySets.Approximations: overapproximate, OctDirections

"""
    function filtered_oscillator(n0::Int=4, time_horizon::Float64=99.,
                                 one_loop_iteration::Bool=false)

Construct the filtered-oscillator model.

### Input

- `n0`                 -- (optional; default: `4`) number of filters
- `time_horizon`       -- (optional; default: `99.`) time horizon
- `one_loop_iteration` -- (optional; default: `false`) option to enforce a
                          single loop (see Notes section below)

### Notes

The option `one_loop_iteration` can be used to add an additional variable `k`.
This variable starts at value `1` and is multiplied by `2` when taking the
transition from location `ℓ₃` to location `ℓ₄`.
In location `ℓ₄` we add the invariant constraint `k ≤ 2.1`.
Thus location `ℓ₄` can only be entered once.
"""
function filtered_oscillator(n0::Int=4,
                             time_horizon::Float64=99.,
                             one_loop_iteration::Bool=false)

    n1 = (one_loop_iteration ? n0 + 1 : n0)
    n = n1 + 2
    z = zeros(n1)

    # transition graph (automaton)
    a = LightAutomaton(4)
    add_transition!(a, 3, 4, 1)
    add_transition!(a, 4, 2, 2)
    add_transition!(a, 2, 1, 3)
    add_transition!(a, 1, 3, 4)

    # common flow
    A = zeros(n, n)
    A[1,1], A[2,2] = -2., -1.
    A[3,1], A[3,3] = 5., -5.
    for i = 4 : n-1
        A[i,i-1], A[i,i] = 5., -5.
    end

    # modes

    # mode 1
    b = [1.4; -0.7; z]
    X = HPolyhedron([HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
                     HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
    m_1 = CACS(A, b, X)

    # mode 2
    b = [-1.4; 0.7; z]
    X = HPolyhedron([HalfSpace([1.0; 0.0; z], 0.0),  # x <= 0
                     HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
    m_2 = CACS(A, b, X)

    # mode 3
    b = [1.4; -0.7; z]
    X = HPolyhedron([HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
                     HalfSpace([-0.714286; -1.0; z], 0.0)])  # 0.714286*x + y >= 0
    m_3 = CACS(A, b, X)

    # mode 4
    b = [-1.4; 0.7; z]
    X = HPolyhedron([HalfSpace([0.714286; 1.0; z], 0.0),  # 0.714286*x + y <= 0
                     HalfSpace([-1.0; 0.0; z], 0.0)])  # x >= 0
    if one_loop_iteration
        # k <= 2 (2.1 to avoid numerical issues)
        addconstraint!(X, HalfSpace([zeros(n-1); 1.], 2.1))
    end
    m_4 = CACS(A, b, X)

    m = [m_1, m_2, m_3, m_4]

    # transitions

    # transition l3 -> l4
    X_l3l4 = HPolyhedron([HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
                          HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
                          HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
    if one_loop_iteration
        A_trans_34 = Matrix(1.0I, n, n)
        A_trans_34[n, n] = 2.  # k' = k * 2
        r1 = ConstrainedLinearMap(A_trans_34, X_l3l4)
    else
        r1 = ConstrainedIdentityMap(n, X_l3l4)
    end

    # transition l4 -> l2
    X_l4l2 = HPolyhedron([HalfSpace([0.714286; 1.0; z], 0.0),  # 0.714286*x + y <= 0
                          HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
                          HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
    r2 = ConstrainedIdentityMap(n, X_l4l2)

    # transition l2 -> l1
    X_l2l1 = HPolyhedron([HalfSpace([1.0; 0.0; z], 0.0),  # x <= 0
                          HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
                          HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
    r3 = ConstrainedIdentityMap(n, X_l2l1)

    # transition l1 -> l3
    X_l1l3 = HPolyhedron([HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
                          HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
                          HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
    r4 = ConstrainedIdentityMap(n, X_l1l3)

    r = [r1, r2, r3, r4]

    # switchings
    s = [HybridSystems.AutonomousSwitching()]

    HS = HybridSystem(a, m, r, s)

    # initial condition in mode 1
    low = [0.2; -0.1; zeros(n1)]
    high = [0.3; 0.1; zeros(n1)]
    if one_loop_iteration
        low[end] = 1.0
        high[end] = 1.0
    end
    X0 = Hyperrectangle(low=low, high=high)

    problem = InitialValueProblem(HS, [(3, X0)])

    options = Options(:T=>time_horizon, :mode=>"reach", :verbosity=>0)

    solver_options = Options(:vars=>1:n, :δ=>0.01, :plot_vars=>[1, 2],
                             :ε_proj=>0.001, :project_reachset=>false)

    return (problem, options, solver_options)
end

"""
    get_projection(sol, projected_dims)

Overapproximate and project a flowpipe.

### Input

- `sol`            -- a flowpipe solution
- `projected_dims` -- an integer vector of length 2 with the variables to project

### Output

A `ReachSolution` containing the overapproximated and projected set.

### Notes

The overapproximation used is octagonal directions.
"""
function get_projection(sol, projected_dims)
    n = sol.options[:n] # system's dimension
    oa = x -> overapproximate(x, OctDirections(n))
    sol_oa = ReachSolution([ReachSet(CartesianProductArray([oa(rs.X)]), rs.t_start, rs.t_end) for rs in sol.Xk], sol.options)
    sol_proj = ReachSolution(project_reach(sol_oa.Xk, projected_dims, n, sol.options), sol.options)
    return sol_proj
end

# single run of filtered oscillator
problem, options, solver_options = filtered_oscillator()
result = solve(problem, options, BFFPSV18(solver_options), ApproximatingDiscretePost());
