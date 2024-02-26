# ===========================================================================
# Spacecraft Rendezvous model
# See https://easychair.org/publications/paper/S2V for a description and
# https://gitlab.com/goranf/ARCH-COMP/blob/master/2018/NLN/C2E2/ARCH18_NLN/rendezvous/c2e2_nonlinear_pass_4d.hyxml
# for a reference model
# ===========================================================================
using SparseArrays, HybridSystems, Reachability, MathematicalSystems,
      MathematicalPredicates, TaylorIntegration

include("eqs_mymul.jl")

function spacecraft_rendezvous(; T=200.0, orderT=10, orderQ=2, abs_tol=1e-10,
                               max_steps=500, plot_vars=[1, 2],
                               project_reachset=false)
    # variables
    x = 1   # x position (negative!)
    y = 2   # y position (negative!)
    vx = 3  # x velocity
    vy = 4  # y velocity
    t = 5   # time

    # number of variables
    n = 4 + 1

    # time of abort
    t_abort = 120.0

    # discrete structure (graph)
    automaton = GraphAutomaton(3)

    # mode 1 ("approaching")
    𝐹 = spacecraft_approaching!
    invariant = HPolyhedron([HalfSpace(sparsevec([x], [1.0], n), -100.0),    # x <= -100
                             HalfSpace(sparsevec([t], [1.0], n), t_abort)])  # t <= t_abort
    m₁ = CBBCS(𝐹, 5, invariant)

    # mode 2 ("rendezvous attempt")
    𝐹 = spacecraft_rendezvous_attempt!
    invariant = HPolyhedron([HalfSpace(sparsevec([x], [-1.0], n), 100.0),           # x >= -100
                             HalfSpace(sparsevec([x], [1.0], n), 100.0),            # x <= 100
                             HalfSpace(sparsevec([y], [-1.0], n), 100.0),           # y >= -100
                             HalfSpace(sparsevec([y], [1.0], n), 100.0),            # y <= 100
                             HalfSpace(sparsevec([x, y], [-1.0, -1.0], n), 141.1),  # x + y >= -141.1
                             HalfSpace(sparsevec([x, y], [1.0, 1.0], n), 141.1),    # x + y <= 141.1
                             HalfSpace(sparsevec([x, y], [1.0, -1.0], n), 141.1),   # -x + y >= -141.1
                             HalfSpace(sparsevec([x, y], [-1.0, 1.0], n), 141.1),   # -x + y <= 141.1
                             HalfSpace(sparsevec([t], [1.0], n), t_abort)])         # t <= t_abort
    m₂ = CBBCS(𝐹, 5, invariant)

    # mode 3 ("aborting")
    𝐹 = spacecraft_aborting!
    invariant = Universe(n)
    m₃ = CBBCS(𝐹, 5, invariant)

    # modes
    modes = [m₁, m₂, m₃]

    # transition 1 -> 2
    add_transition!(automaton, 1, 2, 1)
    guard = HPolyhedron([HalfSpace(sparsevec([x], [-1.0], n), 100.0),           # x >= -100
                         HalfSpace(sparsevec([x], [1.0], n), 100.0),            # x <= 100
                         HalfSpace(sparsevec([y], [-1.0], n), 100.0),           # y >= -100
                         HalfSpace(sparsevec([y], [1.0], n), 100.0),            # y <= 100
                         HalfSpace(sparsevec([x, y], [-1.0, -1.0], n), 141.1),  # x + y >= -141.1
                         HalfSpace(sparsevec([x, y], [1.0, 1.0], n), 141.1),    # x + y <= 141.1
                         HalfSpace(sparsevec([x, y], [1.0, -1.0], n), 141.1),   # -x + y >= -141.1
                         HalfSpace(sparsevec([x, y], [-1.0, 1.0], n), 141.1)])  # -x + y <= 141.1
    t₁ = ConstrainedIdentityMap(n, guard)

    # transition 1 -> 3
    add_transition!(automaton, 1, 3, 2)
    guard = HalfSpace(sparsevec([t], [-1.0], n), -t_abort)  # t >= t_abort
    t₂ = ConstrainedIdentityMap(n, guard)

    # transition 2 -> 3
    add_transition!(automaton, 2, 3, 3)
    guard = HalfSpace(sparsevec([t], [-1.0], n), -t_abort)  # t >= t_abort
    t₃ = ConstrainedIdentityMap(n, guard)

    # transition annotations
    resetmaps = [t₁, t₂, t₃]

    # switching
    switchings = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    X0 = Hyperrectangle([-900.0, -400.0, 0.0, 0.0, 0.0], [25.0, 25.0, 0.0, 0.0, 0.0])
    initial_condition = [(1, X0)]

    𝑃 = InitialValueProblem(ℋ, initial_condition)

    # safety property in "Mode 2"
    velocity = 0.055 * 60.0  # meters per minute
    cx = velocity * cos(π / 8)  # x-coordinate of the octagon's first (ENE) corner
    cy = velocity * sin(π / 8)  # y-coordinate of the octagon's first (ENE) corner

    tan30 = tan(π / 6)
    #=
    octagon = [
        HalfSpace(sparsevec([vx], [1.], n), cx),                 # vx <= cx
        HalfSpace(sparsevec([vx, vy], [1., 1.], n), cy + cx),    # vx + vy <= cy + cx
        HalfSpace(sparsevec([vy], [1.], n), cx),                 # vy <= cx
        HalfSpace(sparsevec([vx, vy], [-1., 1.], n), cy + cx),   # -vx + vy <= cy + cx
        HalfSpace(sparsevec([vx], [-1.], n), cx),                # vx >= -cx
        HalfSpace(sparsevec([vx, vy], [-1., -1.], n), cy + cx),  # -vx - vy <= cy + cx
        HalfSpace(sparsevec([vy], [-1.], n), cx),                # vy >= -cx
        HalfSpace(sparsevec([vx, vy], [1., -1.], n), cy + cx)    # vx - vy <= cy + cx
       ]
    cone = [
        HalfSpace(sparsevec([x], [-1.], n), 100.),          # x >= -100
        HalfSpace(sparsevec([x, y], [tan30, -1.], n), 0.),  # -x tan(30°) + y >= 0
        HalfSpace(sparsevec([x, y], [tan30, 1.], n), 0.),   # -x tan(30°) - y >= 0
       ]
    property_rendezvous = is_contained_in(HPolytope([octagon; cone]))
    =#

    property_rendezvous = (t, x) -> (x[3] <= cx) && (x[3] + x[4] <= cx + cy) && (x[4] <= cx) &&
                                        (-x[3] + x[4] <= cx + cy) && (x[3] >= -cx) &&
                                        (-x[3] - x[4] <= cx + cy) &&
                                        (x[4] >= -cx) && (x[3] - x[4] <= cx + cy) &&
                                        (x[1] >= -100) &&
                                        (-x[1] * tan30 + x[2] >= 0) && (-x[1] * tan30 - x[2] >= 0)

    # safety property in "Passive"
    #=
    target = HPolytope([
        HalfSpace(sparsevec([x], [1.], n), 0.2),   # x <= 0.2
        HalfSpace(sparsevec([x], [-1.], n), 0.2),  # x >= -0.2
        HalfSpace(sparsevec([y], [1.], n), 0.2),   # y <= 0.2
        HalfSpace(sparsevec([y], [-1.], n), 0.2),  # y >= -0.2
       ])
    property_aborting = is_disjoint_from(target)
    =#

    property_aborting = (t, x) -> !((x[1] <= 0.2) && (x[1] >= -0.2) && (x[2] <= 0.2) &&
                                    (x[2] >= -0.2))

    # safety properties
    property = Dict{Int,Function}(1 => (t, x) -> true,
                                  2 => property_rendezvous,
                                  3 => property_aborting)

    # global options
    𝑂 = Options(:T => T, :property => property, :plot_vars => plot_vars,
                :project_reachset => project_reachset, :mode => "reach")

    # algorithm-specific options
    #𝑂jets = Options(:orderT=>orderT, :orderQ=>orderQ, :abs_tol=>abs_tol, :max_steps=>max_steps)

    return 𝑃, 𝑂
end
