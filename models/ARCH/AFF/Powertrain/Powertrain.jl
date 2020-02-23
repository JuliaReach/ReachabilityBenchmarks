using Reachability, MathematicalSystems, MathematicalPredicates, HybridSystems,
      SparseArrays


function print_dynamics(A, b, location_name)
    println("dynamics of location $location_name:")
    for i in 1:size(A, 1)-1  # ignore the last dimension (time)
        print("x_$i' = ")
        for j in 1:size(A, 2)
            if !iszero(A[i,j])
                print("$(A[i,j]) x_$j + ")
            end
        end
        println("$(b[i])\n")
    end
end

"""
    powertrain(θ::Int)

Return the hybrid system that models a mechanical system with backlash
(powertrain) from an automotive problem.

### Input

- `θ`        -- (optional, default: `1`) number of rotating masses
- `X0_scale` -- (optional, default: `1.0`) scaling factor in ``(0, 1]`` for
                modifying the initial states; for value `1.0` we do not modify
                the initial states and let `X0` be just a line segment;
                otherwise we let `X0' = δ ⋅ X0 + ⊕ {(1-δ) · center(X0)}` where
                `δ` is the factor

### Output

Hybrid system representing the powertrain model.

### Notes

The model is based on [1], which is an extension with a parametric number of
rotating masses of the automotive powertrain model presented in [2].
The parameters of the system's additional masses, which can be interpreted as
rotating elements in a gearbox and further powertrain elements, are taken from
[3].

[1] M. Althoff and B. H. Krogh, Avoiding Geometric Intersection Operations in
Reachability Analysis of Hybrid Systems [HSCC '12 Proceedings of the 15th ACM
international conference on Hybrid Systems: Computation and Control Pages
45-54](https://dl.acm.org/citation.cfm?id=2185643).

[2] A. Lagerberg. A benchmark on hybrid control of an automotive powertrain with
backlash. Technical Report R005/2007, Signals and Systems, Chalmers University
of Technology, 2007.

[3] E.-A. M. A. Rabeih. Torsional Vibration Analysis of Automotive Drivelines.
PhD thesis, University of Leeds, 1997.
"""
function powertrain(θ::Int=1; X0_scale::Float64=1.0)
    @assert θ > 0 "θ must be positive, but was $θ"
    @assert (X0_scale > 0.0 && X0_scale <= 1.0) "scale $X0_scale ∉ (0, 1]"

    # activate for printing the dynamics
    display_dynamics = false

    # dimension of state space (last dimension is time)
    n = 2 * θ + 7 + 1

    # constants
    # indices 'm' (ₘ) and 'l' (ₗ) refer to the motor and the load
    # indices 'i' (ᵢ) refer to the numbering of additional rotating masses
    t_init = 0.2  # time to stay in the initial location
    α = 0.03  # backlash size (half of the gap width)
    τ_eng = 0.1  # engine time constant
    γ = 12.0 # gearbox ratio (dimensionless)
    u = 5.0 # requested engine torque
    # moments of inertia [kg m²]
    Jₗ = 140.
    Jₘ = 0.3
    Jᵢ = 100.0  # TODO the paper says 0.01, the SpaceEx model uses 100
    # viscous friction constants [Nm s/rad]
    bₗ = 5.6
    bₘ = 0.0
    bᵢ = 1.0
    # shaft stiffness [Nm/rad]
    kᵢ = 1e5
    kₛ = 1e4
    # PID parameters
    k_P = 0.5  # [Nms/rad]
    k_I = 0.5  # [Nm/rad]
    k_D = 0.5  # [Nms²/rad]

    # additional bloating to get a non-flat guard intersection
    guard_bloating = sqrt(eps(Float64))

    function get_dynamics(kₛ, α, u)
        # physical units of the state variables
        # [x₁] = rad
        # [x₂] = Nm
        # [x₃] = rad
        # [x₄] = rad/s
        # [x₅] = rad
        # [x₆] = rad/s
        # [x₇] = rad/s
        # [x₈] = rad
        # [x₉] = rad/s
        # ...
        # [x_(2θ+6)] = rad
        # [x_(2θ+7)] = rad/s

        # linear dynamics
        A = spzeros(n, n)

        A[1, 7] = 1.0 / γ
        A[1, 9] = -1.0

        A[2, 1] = (-k_I * γ + k_D * kₛ / (γ * Jₘ)) / τ_eng
        A[2, 2] = (-k_D / Jₘ - 1.0) / τ_eng
        A[2, 3] = k_I * γ / τ_eng
        A[2, 4] = k_P * γ / τ_eng
        A[2, 7] = (-k_P + k_D * bₘ /Jₘ) / τ_eng
        A[2, 8] = -k_I * γ / τ_eng

        A[3, 4] = 1.0

        A[5, 6] = 1.0

        A[6, 5] = -kᵢ / Jₗ
        A[6, 6] = -bₗ / Jₗ
        A[6, 2*θ+6] = kᵢ / Jₗ

        A[7, 1] = -kₛ / (Jₘ * γ)
        A[7, 2] = 1.0 / Jₘ
        A[7, 7] = -bₘ / Jₘ

        i = 8
        while i < n-1
            A[i, i+1] = 1.0

            if i == 8
                # x9 has special dynamics
                A[i+1, 1] = kₛ / Jᵢ
                A[i+1, i] = -kᵢ / Jᵢ
            else
                A[i+1, i-2] = kᵢ / Jᵢ
                A[i+1, i] = -2. * kᵢ / Jᵢ
            end
            A[i+1, i+1] = -bᵢ / Jᵢ

            # wrap-around to x5 in the last step
            j = (i == n-2) ? 5 : i+2
            A[i+1, j] = kᵢ / Jᵢ

            i += 2
        end

        # affine vector
        b = spzeros(n)
        b[2] = k_D * (γ * u - kₛ * α / (Jₘ * γ)) / τ_eng
        b[4] = u
        b[7] = kₛ * α / (Jₘ * γ)
        b[9] = -kₛ * α / Jᵢ
        b[n] = 1.0  # time

        return A, b
    end

    # hybrid automaton
    automaton = LightAutomaton(4)

    # negAngle
    A, b = get_dynamics(kₛ, -α, u)
    X = HalfSpace(sparsevec([1], [1.], n), -α)  # x1 <= -α
    m_negAngle = CACS(A, b, X)
    if display_dynamics
        print_dynamics(A, b, "negAngle")
    end

    # deadzone
    A, b = get_dynamics(0., -α, u)
    X = HPolyhedron([HalfSpace(sparsevec([1], [-1.], n), α),  # x1 >= -α
                     HalfSpace(sparsevec([1], [1.], n), α)])  # x1 <= α
    m_deadzone = CACS(A, b, X)
    if display_dynamics
        print_dynamics(A, b, "deadzone")
    end

    # posAngle
    A, b = get_dynamics(kₛ, α, u)
    X = HalfSpace(sparsevec([1], [-1.], n), -α + guard_bloating)  # x1 >= α - ε
    m_posAngle = CACS(A, b, X)
    if display_dynamics
        print_dynamics(A, b, "posAngle")
    end

    # negAngleInit
    A, b = get_dynamics(kₛ, -α, -u)
    X = HalfSpace(sparsevec([n], [1.], n), t_init)  # t <= t_init
    m_negAngleInit = CACS(A, b, X)
    if display_dynamics
        print_dynamics(A, b, "negAngleInit")
    end

    # modes
    modes = [m_negAngle, m_deadzone, m_posAngle, m_negAngleInit]

    # transition negAngleInit -> negAngle
    add_transition!(automaton, 4, 1, 1)
    guard = HalfSpace(sparsevec([n], [-1.], n), -t_init + guard_bloating)  # t >= t_init - ε
    r_41 = ConstrainedIdentityMap(n, guard)

    # transition negAngle -> deadzone
    add_transition!(automaton, 1, 2, 2)
    guard = HalfSpace(sparsevec([1], [-1.], n), α)  # x1 >= -α
    r_12 = ConstrainedIdentityMap(n, guard)

    # transition deadzone -> posAngle
    add_transition!(automaton, 2, 3, 3)
    guard = HalfSpace(sparsevec([1], [-1.], n), -α)  # x1 >= α
    r_23 = ConstrainedIdentityMap(n, guard)

    # TODO the SpaceEx model does not contain the following transitions
#     # transition deadzone -> negAngle
#     add_transition!(automaton, 2, 1, 4)
#     guard = HalfSpace(sparsevec([1], [1.], n), -α)  # x1 <= -α
#     r_21 = ConstrainedIdentityMap(n, guard)
#     # transition posAngle -> deadzone
#     add_transition!(automaton, 3, 2, 5)
#     guard = HalfSpace(sparsevec([1], [1.], n), α)  # x1 <= α
#     r_32 = ConstrainedIdentityMap(n, guard)

    # transition annotations
    resetmaps = [r_41, r_12, r_23]

    # switching
    switchings = [HybridSystems.AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition in mode 1
    c = Vector{Float64}(undef, n)
    g = Vector{Float64}(undef, n)
    c[1:7] = [-0.0432, -11., 0., 30., 0., 30., 360.]
    g[1:7] = [0.0056, 4.67, 0., 10., 0., 10., 120.]
    i = 8
    while i < n
        c[i] = -0.0013
        g[i] = 0.0006
        i += 1
        c[i] = 30.
        g[i] = 10.
        i += 1
    end
    c[n] = 0.0
    g[n] = 0.0
    if X0_scale < 1.0
        g = X0_scale * g
    end
    X0 = Zonotope(c, hcat(g))
    initial_condition = [(4, X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    # safety property
    property_2 = is_disjoint_from(HalfSpace(sparsevec([1], [1.], n), -α))  # x1 <= -α
    property_3 = is_disjoint_from(HalfSpace(sparsevec([1], [1.], n), α))   # x1 <= α
    property = Dict(2 => property_2, 3 => property_3)

    # default options
    options = Options(:T=>2.0, :property=>property)

    return (system, options)
end

function run_powertrain(system, options)
    opC = BFFPSV18(:δ => 0.0005, :assume_sparse => true)
    opD = LazyDiscretePost(:lazy_R⋂I => true, :lazy_R⋂G => true)
    options[:mode] = "check"
    options[:plot_vars] = [1, 3]

    solve(system, options, opC, opD)
end
