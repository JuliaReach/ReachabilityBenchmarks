# ===========================================================================
# Spacecraft Rendezvous model
# See https://easychair.org/publications/paper/S2V for a description and
# https://gitlab.com/goranf/ARCH-COMP/blob/master/2018/NLN/C2E2/ARCH18_NLN/rendezvous/c2e2_nonlinear_pass_4d.hyxml
# for a reference model
# ===========================================================================
using SparseArrays, HybridSystems, MathematicalSystems, LazySets, Reachability, TaylorIntegration
using Reachability: solve

const μ = 3.986e14 * 60^2
const r = 42164e3
const mc = 500.0
const n = sqrt(μ / r^3)

# dynamics
@taylorize function spacecraft_approaching!(t, x, dx, K)
    local rc = sqrt((r + x[1])^2 + x[2]^2)
    local uxy = K * x

    # x' = vx
    dx[1] = x[3]

    # y' = vy
    dx[2] = x[4]

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    dx[3] = (n^2*x[1] + 2*(n*x[4])) + ((μ/(r^2) - μ/(rc^3)*(r + x[1])) + uxy[1]/mc)

    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    dx[4] = (n^2*x[2] - 2*(n*x[3])) - (μ/(rc^3)*x[2] - uxy[2]/mc)

    return dx
end

# dynamics in the aborting mode
@taylorize function spacecraft_aborting!(t, x, dx)
    local rc = sqrt((r + x[1])^2 + x[2]^2)

    # x' = vx
    dx[1] = x[3]

    # y' = vy
    dx[2] = x[4]

    # vx' = n²x + 2n*vy + μ/(r^2) + μ/(rc^3)*(r+x)
    dx[3] = (n^2*x[1] + 2*n*x[4]) + μ/(r^2) + μ/(rc^3)*(r + x[1])

    # vy' = n²y - 2n*vx - μ/(rc^3)y
    dx[4] = (n^2*x[2] - 2*n*x[3]) - μ/(rc^3)*x[2]

    return dx
end

function spacecraft_rendezvous()
    # number of variables
    n = 4

    # time of abort
    t_abort = 120.0

    # discrete structure (graph)
    automaton = LightAutomaton(3)

    # mode 1 ("approaching")
    K₁ = [-28.8287 0.1005 -1449.9754 0.0046;
          -0.087 -33.2562 0.00462 -1451.5013]
    𝐹 = (t, x, dx) -> spacecraft_approaching!(t, x, dx, K₁)
    invariant = HalfSpace(sparsevec([x], [1.], n), -100.)  # x <= -100
    m₁ = CBBCS(𝐹, 4, invariant)

    # mode 2 ("rendezvous attempt")
    K₂ = [-288.0288 0.1312 -9614.9898 0.0;
          -0.1312 -288.0 0.0 -9614.9883]
    𝐹 = (t, x, dx) -> spacecraft_approaching!(t, x, dx, K₂)
    invariant = HPolyhedron([
        HalfSpace(sparsevec([x], [-1.], n), 100.),           # x >= -100
        HalfSpace(sparsevec([x], [1.], n), 100.),            # x <= 100
        HalfSpace(sparsevec([y], [-1.], n), 100.),           # y >= -100
        HalfSpace(sparsevec([y], [1.], n), 100.),            # y <= 100
        HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1),  # x + y >= -141.1
        HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1),    # x + y <= 141.1
        HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1),   # -x + y >= -141.1
        HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1)    # -x + y <= 141.1
       ])
    m₂ = CBBCS(𝐹, 4, invariant)

    # mode 3 ("aborting")
    𝐹 = spacecraft_aborting!
    invariant = Universe(n)
    m₃ = CBBCS(𝐹, 4, invariant)

    # modes
    modes = aborting ? [m₁, m₂, m₃]

    # transition 1 -> 2
    add_transition!(automaton, 1, 2, 1)
    guard = HPolyhedron([
        HalfSpace(sparsevec([y], [-1.], n), 100.),           # y >= -100
        HalfSpace(sparsevec([x, y], [-1., -1.], n), 141.1),  # x + y >= -141.1
        HalfSpace(sparsevec([x], [-1.], n), 100.),           # x >= -100
        HalfSpace(sparsevec([x, y], [-1., 1.], n), 141.1),   # -x + y <= 141.1
        HalfSpace(sparsevec([y], [1.], n), 100.),            # y <= 100
        HalfSpace(sparsevec([x, y], [1., 1.], n), 141.1),    # x + y <= 141.1
        HalfSpace(sparsevec([x], [1.], n), 100.),            # x <= 100
        HalfSpace(sparsevec([x, y], [1., -1.], n), 141.1),   # -x + y >= -141.1
       ])
    t₁ = ConstrainedIdentityMap(n, guard)

    # transition 1 -> 3
    add_transition!(automaton, 1, 3, 2)
    guard = HalfSpace(sparsevec([t], [-1.], n), -t_abort)  # t >= t_abort
    t₂ = ConstrainedIdentityMap(n, guard)

    # transition 2 -> 3
    add_transition!(automaton, 2, 3, 3)
    guard = HalfSpace(sparsevec([t], [-1.], n), -t_abort)  # t >= t_abort
    t₃ = ConstrainedIdentityMap(n, guard)

    # transition annotations
    resetmaps = [t₁, t₂, t₃]

    # switching
    switchings = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    X0 = Hyperrectangle([-900., -400., 0., 0.], [25., 25., 0., 0.])
    initial_condition = [(1, X0)]

    𝑃 = InitialValueProblem(ℋ, initial_condition)

    # TODO
    property = (t, x) -> true

    𝑂 = Options(:T=>200.0, :property=>property)

    return 𝑃, 𝑂
end

function spacecraft_TMJets(property; T=200.0, orderT=10, orderQ=2, abs_tol=1e-10, max_steps=500)
    𝑃, 𝑂 = spacecraft_rendezvous()
    𝑂jets = Options(:orderT=>orderT, :orderQ=>orderQ, :abs_tol=>abs_tol, :max_steps=>max_steps)
    return 𝑃, 𝑂, 𝑂jets
end
