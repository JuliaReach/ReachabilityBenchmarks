# =================================================
# Platooning model
# See https://easychair.org/publications/paper/3QLs
# =================================================
using LinearAlgebra, SparseArrays, HybridSystems, MathematicalSystems, LazySets,
      Reachability
using LazySets: HalfSpace  # resolve name-space conflicts with Polyhedra

"""
    platooning(;
               [deterministic_switching]::Bool=true,
               [time_horizon]::Float64=20.,
               [allowed_distance]::Float64=50.)

Construct the platooning model consisting of three vehicles.

### Input

- `deterministic_switching` -- (optional, default: `true`) flag to use
                               deterministic transitions
- `time_horizon`            -- (optional, default: `20.`) time horizon
- `allowed_distance`        -- (optional, default: `50.`) specification of
                               allowed distance between the vehicles
"""
function platooning(;
                    deterministic_switching::Bool=true,
                    time_horizon::Float64=20.,
                    allowed_distance::Float64=50.)
    # three variables for each vehicle, (ei, d(et)/dt, ai) for
    # (spacing error, relative velocity, speed), and the last dimension is time
    n = 9 + 1

    # constants
    c1 = c2 = 5.   # clock constraints
    tb = 10.       # lower bound for loss of communication
    tc = tr = 20.  # upper bound for loss of communication (tc) and reset time (tr)

    # additional bloating to get a non-flat guard intersection
    guard_bloating = sqrt(eps(Float64))

    # transition graph
    automaton = LightAutomaton(2)
    add_transition!(automaton, 1, 2, 1)
    add_transition!(automaton, 2, 1, 2)

    # common inputs
    B = sparse([2, n], [1, 2], [1., 1.], n, 2)
    U = Hyperrectangle(low=[-9., 1.], high=[1., 1.])  # acceleration of the lead vehicle + time

    # mode 1 ("connected")
    A = zeros(n, n)
    for (i, j) in [(1, 2), (4, 5), (5, 3), (7, 8), (8, 6)]
        A[i, j] = 1.
    end
    for (i, j) in [(2, 3), (5, 6), (8, 9)]
        A[i, j] = -1.
    end
    A[3, :] = [1.6050, 4.8680, -3.5754, -0.8198, 0.4270, -0.0450, -0.1942,  0.3626, -0.0946, 0.]
    A[6, :] = [0.8718, 3.8140, -0.0754,  1.1936, 3.6258, -3.2396, -0.5950,  0.1294, -0.0796, 0.]
    A[9, :] = [0.7132, 3.5730, -0.0964,  0.8472, 3.2568, -0.0876,  1.2726,  3.0720, -3.1356, 0.]
    invariant = deterministic_switching ?
        HalfSpace(sparsevec([n], [1.], n), c1) :
        Universe(n)
    m_1 = ConstrainedLinearControlContinuousSystem(A, B, invariant, U)

    # mode 2 ("not connected/connection broken")
    A = copy(A)
    A[3, 4:9] = zeros(6)
    A[6, 1:3] = A[6, 7:9] = zeros(3)
    invariant = deterministic_switching ?
        HalfSpace(sparsevec([n], [1.], n), c2) :
        Universe(n)
    m_2 = ConstrainedLinearControlContinuousSystem(A, B, invariant, U)

    # modes
    modes = [m_1, m_2]

    # common reset
    reset = Dict(n => 0.)

    # transition l1 -> l2
    # (using a hyperplane in the deterministic case causes floating-point issues)
    guard = deterministic_switching ?
        HalfSpace(sparsevec([n], [-1.], n), -c1 + guard_bloating) :
        HPolyhedron([HalfSpace(sparsevec([n], [-1.], n), -tb),
                     HalfSpace(sparsevec([n], [1.], n), tc)])
    t1 = ConstrainedResetMap(n, guard, reset)

    # transition l2 -> l1
    guard = deterministic_switching ?
        HalfSpace(sparsevec([n], [-1.], n), -c2 + guard_bloating) :
        HalfSpace(sparsevec([n], [1.], n), tr)
    t2 = ConstrainedResetMap(n, guard, reset)

    # transition annotations
    resetmaps = [t1, t2]

    # switching
    switchings = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition in mode 1
    X0 = Singleton(zeros(n))
    initial_condition = [(1, X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    # safety property
    d1 = zeros(n); d1[1] = -1.  # x1 >= -dmin
    d4 = zeros(n); d4[4] = -1.  # x4 >= -dmin
    d7 = zeros(n); d7[7] = -1.  # x7 >= -dmin
    property = Conjunction(
        [SafeStatesProperty(HalfSpace(d, allowed_distance)) for d in [d1, d4, d7]])

    # default options
    options = Options(:T=>time_horizon, :property=>property)

    return (system, options)
end

function run_platooning(system, options)
    opC = BFFPSV18(:δ => 0.01, :assume_sparse => true)
    opD = LazyDiscretePost(:lazy_R⋂I => true, :lazy_R⋂G => false)
    options[:mode] = "check"
    options[:plot_vars] = [0, 1]

    solve(system, options, opC, opD)
end


function run_platooning(system, options, opD, opC)
    options[:mode] = "check"
    println("upd")
    return solve(system, options, opC, opD)
end
