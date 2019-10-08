# ===========================================================================
# Spacecraft Rendezvous model
# See https://easychair.org/publications/paper/S2V for a description and
# https://gitlab.com/goranf/ARCH-COMP/tree/master/2018/AFF/SpaceEx/Rendezvous
# for a reference model
# ===========================================================================
using SparseArrays, HybridSystems, MathematicalSystems, LazySets, Reachability
using LazySets: HalfSpace  # resolve name-space conflicts with Polyhedra

"""
    spacecraft(; abort_time::Union{Float64, Vector{Float64}}=-1.)

Construct the spacecraft-rendezvous model.

### Input

- `abort_time` -- (optional, default: `-1.`) time of abort; can be either a
                  number (time point) or an interval (time range); an interval
                  is given as a vector `[l, u]` containing the lower (`l`) and
                  the upper (`u`) bound; a negative number is interpreted as a
                  *no-abort* scenario
"""
function spacecraft(; abort_time::Union{Float64, Vector{Float64}}=120.)
    # variables
    x = 1  # x position (negative!)
    y = 2  # y position (negative!)
    vx = 3  # x velocity
    vy = 4  # y velocity
    t = 5   # time

    # number of variables
    n = 4 + 1
    ε = 1e-6
    # flag for activating the "rendezvous abort" scenario
    aborting = true
    if abort_time isa Number
        if abort_time < 0.
            aborting = false
        else
            t_abort_lower = abort_time
            t_abort_upper = abort_time
        end
    else
        @assert length(abort_time) == 2 "abort time must be a point or an interval"
        t_abort_lower = abort_time[1]
        t_abort_upper = abort_time[2]
    end

    # discrete structure (graph)
    automaton = LightAutomaton(aborting ? 3 : 2)

    # common vector of affine dynamics
    b = sparsevec([t], [1.], n)

    # mode 1 ("approaching")
    A = spzeros(n, n)
    A[x, vx] = 1.
    A[y, vy] = 1.
    A[vx, x] = -0.057599765881773
    A[vx, y] = 0.000200959896519766
    A[vx, vx] = -2.89995083970656
    A[vx, vy] = 0.00877200894463775
    A[vy, x] = -0.000174031357370456
    A[vy, y] = -0.0665123984901026
    A[vy, vx] = -0.00875351105536225
    A[vy, vy] = -2.90300269286856
    invariant = HalfSpace(sparsevec([x], [1.], n), -100.)  # x <= -100
    if aborting
        invariant = HPolyhedron([
            invariant,
            HalfSpace(sparsevec([t], [1.], n), t_abort_upper)  # t <= t_abort_upper
           ])
    end
    m_1 = CACS(A, b, invariant)

    # mode 2 ("rendezvous attempt")
    A = spzeros(n, n)
    A[x, vx] = 1.
    A[y, vy] = 1.
    A[vx, x] = -0.575999943070835
    A[vx, y] = 0.000262486079431672
    A[vx, vx] = -19.2299795908647
    A[vx, vy] = 0.00876275931760007
    A[vy, x] = -0.000262486080737868
    A[vy, y] = -0.575999940191886
    A[vy, vx] = -0.00876276068239993
    A[vy, vy] = -19.2299765959399
    invariant = HPolyhedron([
        HalfSpace(sparsevec([x], [-1.], n), 100. + ε),           # x >= -100
        HalfSpace(sparsevec([x], [1.], n), 100. + ε),            # x <= 100
        HalfSpace(sparsevec([y], [-1.], n), 100. + ε),           # y >= -100
        HalfSpace(sparsevec([y], [1.], n), 100. + ε),            # y <= 100
        HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1+ ε),  # x + y >= -141.1
        HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1+ ε),    # x + y <= 141.1
        HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1+ ε),   # -x + y >= -141.1
        HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1+ ε)    # -x + y <= 141.1
       ])
    if aborting
        addconstraint!(invariant, HalfSpace(sparsevec([t], [1.], n), t_abort_upper+ ε))  # t <= t_abort_upper
    end
    m_2 = CACS(A, b, invariant)

    # mode 3 ("aborting")
    A = spzeros(n, n)
    A[x, vx] = 1.
    A[y, vy] = 1.
    A[vx, x] = 0.0000575894721132
    A[vx, vy] = 0.00876276
    A[vy, vx] = -0.00876276
    invariant = Universe(n)
    m_3 = CACS(A, b, invariant)

    # modes
    modes = aborting ? [m_1, m_2, m_3] : [m_1, m_2]

    # transition 1 -> 2
    add_transition!(automaton, 1, 2, 1)
    guard = HPolyhedron([
        HalfSpace(sparsevec([x], [-1.], n), 100. + ε),           # x >= -100
        HalfSpace(sparsevec([x], [1.], n), 100. + ε),            # x <= 100
        HalfSpace(sparsevec([y], [-1.], n), 100. + ε),           # y >= -100
        HalfSpace(sparsevec([y], [1.], n), 100. + ε),            # y <= 100
        HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1+ ε),  # x + y >= -141.1
        HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1+ ε),    # x + y <= 141.1
        HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1+ ε),   # -x + y >= -141.1
        HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1+ ε)    # -x + y <= 141.1
       ])
    t1 = ConstrainedIdentityMap(n, guard)

    # TODO the SpaceEx model is missing the transition from 2 to 1

    if aborting
        add_transition!(automaton, 1, 3, 2)
        guard = HPolyhedron([HalfSpace(sparsevec([t], [-1.], n), -t_abort_lower+ ε)])  # t >= t_abort_lower
        t2 = ConstrainedIdentityMap(n, guard)

        add_transition!(automaton, 2, 3, 3)
        guard = HPolyhedron([HalfSpace(sparsevec([t], [-1.], n), -t_abort_lower+ ε)])  # t >= t_abort_lower
        t3 = ConstrainedIdentityMap(n, guard)
    end

    # transition annotations
    resetmaps = aborting ? [t1, t2, t3] : [t1]

    # switching
    switchings = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition in mode 1
    X0 = Hyperrectangle([-900., -400., 0., 0., 0.],
                        [25., 25., 0., 0., 0.])
    initial_condition = [(1, X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    # safety property in "Mode 2"
    velocity = 0.055 * 60.  # meters per minute
    cx = velocity * cos(π / 8)  # x-coordinate of the octagon's first (ENE) corner
    cy = velocity * sin(π / 8)  # y-coordinate of the octagon's first (ENE) corner
    octagon = [
        HalfSpace(sparsevec([vx], [-1.], n), cx),                # vx >= -cx
        HalfSpace(sparsevec([vx], [1.], n), cx),                 # vx <= cx
        HalfSpace(sparsevec([vy], [-1.], n), cx),                # vy >= -cx
        HalfSpace(sparsevec([vy], [1.], n), cx),                 # vy <= cx
        HalfSpace(sparsevec([vx, vy], [1., 1.], n), cy + cx),    # vx + vy <= cy + cx
        HalfSpace(sparsevec([vx, vy], [1., -1.], n), cy + cx),   # vx - vy <= cy + cx
        HalfSpace(sparsevec([vx, vy], [-1., 1.], n), cy + cx),   # -vx + vy <= cy + cx
        HalfSpace(sparsevec([vx, vy], [-1., -1.], n), cy + cx)   # -vx - vy <= cy + cx
       ]
    tan30 = tan(π/6)
    cone = [
        HalfSpace(sparsevec([x], [-1.], n), 100.),          # x >= -100
        HalfSpace(sparsevec([x, y], [tan30, -1.], n), 0.),  # -x tan(30°) + y >= 0
        HalfSpace(sparsevec([x, y], [tan30, 1.], n), 0.),   # -x tan(30°) - y >= 0
       ]
    property_rendezvous = SafeStatesProperty(HPolyhedron([octagon; cone]))
    # safety property in "aborting"
    target = HPolyhedron([
        HalfSpace(sparsevec([x], [-1.], n), 0.2),  # x >= -0.2
        HalfSpace(sparsevec([x], [1.], n), 0.2),   # x <= 0.2
        HalfSpace(sparsevec([y], [-1.], n), 0.2),  # y >= -0.2
        HalfSpace(sparsevec([y], [1.], n), 0.2),   # y <= 0.2
       ])
    property_aborting = BadStatesProperty(target)
    # safety properties
    property = Dict{Int, Property}(2 => property_rendezvous,
                                   3 => property_aborting)

    # default options
    options = Options(:T=>240., :property=>property)

    return (system, options)
end

function run_spacecraft(system, options, opC, opD)
    options[:mode] = "check"

    solve(system, options, opC, opD)
end
