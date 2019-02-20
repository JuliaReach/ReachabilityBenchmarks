# ================================================
# Gearbox model
# See https://easychair.org/publications/paper/cwl
# ================================================
using LinearAlgebra, SparseArrays, HybridSystems, MathematicalSystems, LazySets,
      Reachability
using LazySets: HalfSpace  # resolve name-space conflicts with Polyhedra

"""
    gearbox()

Construct the gearbox model.
"""
function gearbox()
    # variables
    vx = 1  # x velocity
    vy = 2  # y velocity
    px = 3  # x position
    py = 4  # y position
    𝐼 = 5  # accumulated impulse (italic version to avoid a clash with 'I')
    t = 6   # time

    # number of variables
    n = 5 + 1

    # constants
    Fs = 70.0  # shifting force on the sleeve [N]
    Tf = 1.0  # resisting moments on the second gear [N⋅m]
    # TODO where do the following values come from?
    ms = 3.2  # mass of the sleeve
    Rs = 0.08  # radius of the sleeve
    Jg₂ = 0.7  # inertia of the second gear
    Δp = -0.003  # horizontal (px) distance
    θ = 0.628318530717959  # included angle of the second gear [°]
    mg₂ = 18.1  # mass of the second gear
    ζ = 0.9  # coefficient of restitution

    # additional bloating to get a non-flat guard intersection
    guard_bloating = sqrt(eps(Float64))

    # discrete structure (graph)
    automaton = LightAutomaton(2)
    add_transition!(automaton, 1, 1, 1)
    add_transition!(automaton, 1, 1, 2)
    add_transition!(automaton, 1, 2, 3)

    # common inputs
    BB = Diagonal(ones(n))
    U = Singleton([1.])

    # mode 1 ("free")
    A = zeros(n, n)
    B = zeros(n, 1)
    A[px, vx] = 1.
    A[py, vy] = 1.
    B[vx, 1] = Fs / ms
    B[vy, 1] = - (Rs * Tf) / Jg₂
    B[t, 1] = 1.
    invariant = HalfSpace(sparsevec([px], [1.], n), Δp)
    # TODO The SpaceEx model adds more constraints, possibly to help with the
    # guard intersection:
    # py <= -px * tan(θ)  &  py >= px * tan(θ)
    m_1 = ConstrainedLinearControlContinuousSystem(A, BB, invariant, B*U)

    # mode 2 ("meshed")
    A = zeros(n, n)
    B = zeros(n, 1)
    invariant = Universe(n)
    m_2 = ConstrainedLinearControlContinuousSystem(A, BB, invariant, B*U)

    # modes
    modes = [m_1, m_2]

    # common assignment matrix (requires individual modifications)
    A_template = zeros(n, n)
    for i in 3:6
        A_template[i, i] = 1.
    end
    denominator = ms * cos(θ)^2 + mg₂ * sin(θ)^2
    A_template[vx, vx] = (ms * cos(θ)^2 - mg₂ * ζ * sin(θ)^2) / denominator
    A_template[vx, vy] = (-(ζ + 1.) * mg₂ * sin(θ) * cos(θ)) / denominator
    A_template[vy, vx] = (-(ζ + 1.) * ms * sin(θ) * cos(θ)) / denominator
    A_template[vy, vy] = (mg₂ * sin(θ)^2 - ms * ζ * cos(θ)^2) / denominator
    A_template[𝐼, vx] = ((ζ + 1.) * ms * mg₂ * sin(θ)) / denominator
    A_template[𝐼, vy] = ((ζ + 1.) * ms * mg₂ * cos(θ)) / denominator

    # transition l1 -> l1
    # TODO what happened to the term '2nb' and the whole second constraint in the paper?
    guard = HPolyhedron([
        HalfSpace(sparsevec([px, py], [-tan(θ), -1.], n), 0.),  # px * tan(θ) + py >= 0
        HalfSpace(sparsevec([vx, vy], [-sin(θ), -cos(θ)], n), 0.)  # vx * sin(θ) + vy * cos(θ) >= 0
        ])
    A = copy(A_template)
    t1 = ConstrainedLinearMap(A, guard)

    # transition l1 -> l1
    # TODO same remark as with the other guard
    guard = HPolyhedron([
        HalfSpace(sparsevec([px, py], [-tan(θ), 1.], n), 0.),  # -px * tan(θ) + py <= 0
        HalfSpace(sparsevec([vx, vy], [-sin(θ), cos(θ)], n), 0.)  # vx * sin(θ) - vy * cos(θ) >= 0
        ])
    A = copy(A_template)
    A[vx, vy] *= -1.
    A[vy, vx] *= -1.
    A[𝐼, vy] *= -1.
    t2 = ConstrainedLinearMap(A, guard)

    # transition l1 -> l2
    guard = HalfSpace(sparsevec([px], [-1.], n), -Δp + guard_bloating)  # px >= Δp
    A = copy(A_template)
    A[vx, vx] = 0.
    A[vx, vy] = 0.
    A[vy, vx] = 0.
    A[vy, vy] = 0.
    A[𝐼, vx] = A[𝐼, vy] = ms
    t3 = ConstrainedLinearMap(A, guard)
    # TODO The SpaceEx model uses four transitions for this one, possibly to
    # help with the signs:
    # * guard = px >= Δp & vx >= 0 & vy >= 0
    #   assignment = I:=I+ms*vx+ms*vy & vx:=0 & vy:=0
    # * guard = px >= Δp & vx >= 0 & vy <= 0
    #   assignment = I:=I+ms*vx-ms*vy & vx:=0 & vy:=0
    # * guard = px >= Δp & vx <= 0 & vy >= 0
    #   assignment = I:=I-ms*vx+ms*vy & vx:=0 & vy:=0
    # * guard = px >= Δp & vx <= 0 & vy <= 0
    #   assignment = I:=I-ms*vx-ms*vy & vx:=0 & vy:=0

    # transition annotations
    resetmaps = [t1, t2, t3]

    # switching
    switchings = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition in mode 1
    X0 = Hyperrectangle([0., 0., -0.0167, 0.003,  0., 0.],
                        [0., 0.,  0.0001, 0.0001, 0., 0.])
    initial_condition = [(1, X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    # safety property
    cond_free = SafeStatesProperty(HalfSpace([0., 0., 0., 0., 0., 1.], 0.2))    # t <= 0.2
    cond_global = SafeStatesProperty(HalfSpace([0., 0., 0., 0., 1., 0.], 20.))  # I <= 20
    property_free = Conjunction([cond_free, cond_global])
    property_meshed = cond_global
    property = Dict(1 => property_free, 2 => property_meshed)

    # default options
    options = Options(:T=>0.5, :property=>property)

    return (system, options)
end

function run_gearbox(system, options)
    opC = BFFPSV18(:δ => 0.001)
    opD = LazyDiscretePost(:lazy_R⋂I => true, :lazy_R⋂G => false)
    options[:verbosity] = "info"
    options[:mode] = "check"
    options[:plot_vars] = [3, 4]
    options[:project_reachset] = true

    solve(system, options, opC, opD)
end
