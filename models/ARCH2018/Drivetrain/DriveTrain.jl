using MathematicalSystems, HybridSystems, LazySets

# compatibility of julia versions
if VERSION >= v"0.7"
    using LinearAlgebra, SparseArrays
end

"""
    drivetrain(θ::Int)::HybridSystem

Return the hybrid system that models a mechanical system with backlash
(powertrain) from an automotive drivetrain problem.

### Input

- `θ` -- (optional, default: `1`) number of rotating masses

### Output

Hybrid system representing the powertrain model.

### Notes

The model is based on [1], which is an extension with a parametric number of
rotationg massesof the autmotive powertrain model presented in [2].
The parameters of the system's additional masses, which can be interpreted as
rotating elements in a gearbox and further drivetrain elements, are taken from
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
function drivetrain(θ::Int=1)::HybridSystem

    # dimension of state space
    n = 2 * θ + 7

    # =========
    # constants
    # =========
    α = 0.03 # backlash size (half gap width)
    α_neg = -0.03
    τ_eng = 0.1 # engine time constant
    γ = 12. # the gearbox ratio
    u = 5 # requested engine torque

    # indices m and l refer to the motor and the load
    # numbered indices refer to the numbering of additional rotating masses,
    # which are currently generalized by i

    # viscous friction constants by b [Nm/(s/rad)]
    b_l = 5.6
    b_m = 0.
    b_i = 1.
    # moments of inertia are denoted by J [kg m²]
    J_l = 140.
    J_m = 0.3
    J_i = 0.01

    # shaft stiffness by k [Nm/rad]
    k_i = 100000.
    k_s = 10000.
    k_s_zero = 0.

    # PID parameters
    k_P = 0.5
    k_I = 0.5
    k_D = 0.5

    # =========
    # dynamics
    # =========
    function get_dynamics(k_s, α)
        # common flow
        A = zeros(n, n) # use spzeros? => see #40 in MathematicalSystems.jl

        J_arr = fill(J_i, θ)
        b_arr = fill(b_i, θ)
        k_arr = fill(k_i, θ)

        A[1,7] = 1.0/γ
        A[1,9] =  -1.
        A[2,1] = -(1.0/τ_eng)*k_I*γ+k_D*(1.0/(J_m*γ*τ_eng)*k_s)
        A[2,2] = -k_D*(1.0/J_m - 1)/τ_eng
        A[2,3] = (1.0/τ_eng)*k_I*γ
        A[2,4] = (1.0/τ_eng)*k_P*γ
        A[2,7] = (-k_P + k_D*(1.0/J_m)*b_m)/τ_eng
        A[2,8] = -(1.0/τ_eng)*k_I

        A[3,4] = 1.

        A[5,6] = 1.

        A[6,5] = -(1.0/J_l) * k_arr[θ]
        A[6,6] =  -(1.0/J_l) * b_l
        A[6,2 * θ+6] = (1.0/J_l) * k_arr[θ]

        A[7,1] = -(1.0/(J_m*γ))*k_s
        A[7,2] = 1.0/J_m
        A[7,7] = -(1.0/J_m)*b_m

        i = 10
        if (θ >= 1)
            A[8, 9] = 1.
            A[9,1] = (1.0/J_arr[1])*k_s
            A[9,8] = -(1.0/J_arr[1])*k_arr[1]
            A[9,9] = -(1.0/J_arr[1])*b_arr[1]
            if (ϴ > 1)
                A[9,10] = (1.0/J_arr[1])*k_arr[1]
            else
                A[9,5] = (1.0/J_arr[1])*k_arr[1]
            end
        end
        el = 2 # returns the index of additional rotating mass
        while i < n
            A[i,i+1] = 1.
            A[i+1,2*el+4] = (1.0/J_arr[el])*k_arr[el-1]
            A[i+1,2*el+6] = -(1.0/J_arr[el])*(k_arr[el-1]+k_arr[el])
            A[i+1,2*el+7] =  -(1.0/J_arr[el])*b_arr[el]
            if el == ϴ
                A[i+1,5] = (1.0/J_arr[el])*k_arr[el]
            else
                A[i+1,i+2] = (1.0/J_arr[el])*k_arr[el]
            end
            el +=1
            i += 2
        end
        return A
    end

    function get_b(k_s, α)
        x2_u = (k_D*γ*u-k_D*
                (1.0/(J_m*γ))*k_s*α)/τ_eng
        x4_u = u
        x7_u = (1.0/(J_m*γ))*k_s*α
        x9_u = -(1.0/J_i)*k_s*α

        return sparsevec([2, 4, 7, 9], [x2_u, x4_u, x7_u, x9_u], n)
    end

    # =============================
    # transition graph (automaton)
    # =============================
    automaton = LightAutomaton(3);

    add_transition!(automaton, 1, 2, 1);
    add_transition!(automaton, 2, 1, 1);
    add_transition!(automaton, 2, 3, 1);
    add_transition!(automaton, 3, 2, 1);

    # common U
    U = Singleton([1.0])

    # common resets
    A_trans = Matrix{Float64}(I, n, n) # use UniformScaling n * I ? not accepted since it doesn't subtype <:AbstractArray{T,2}

    z = zeros(n)

    # negAngle
    A = get_dynamics(k_s, α_neg)
    B = get_b(k_s, α_neg)

    # identity matrix
    E = Matrix{Float64}(I, n, n)

    X = HPolyhedron([HalfSpace([1.; z], α_neg)]) # x <= -α
    m_negAngle = ConstrainedLinearControlContinuousSystem(A, E, X, B*U)

    # transition negAngle -> deadzone
    X_l1l2 = HPolyhedron([HalfSpace([-1.; z], α)])  # x >= -0.03
    r1 = ConstrainedLinearDiscreteSystem(A_trans, X_l1l2);

    # deadzone
    A = get_dynamics(k_s_zero, α_neg)
    B = get_b(k_s_zero, α_neg)
    X = HPolyhedron([HalfSpace([-1.; z], α),  # x >= -α
                     HalfSpace([1.; z], α)])  # x <= 0.03
    m_deadzone = ConstrainedLinearControlContinuousSystem(A, E, X, B*U)

    # transition deadzone -> negAngle
    X_l2l1 = HPolyhedron([HalfSpace([1.; z], α_neg)])  # x <= -0.03
    r2 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l1)

    # transition deadzone -> posAngle
    X_l2l3 = HPolyhedron([HalfSpace([-1.; z], α_neg)])  # x >= 0.03
    r3 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l3)

    #posAngle
    A = get_dynamics(k_s, α)
    B = get_b(k_s, α)
    X = HPolyhedron([HalfSpace([-1.; z], α_neg)])  # x >= α (2.1 for numerical issues)
    m_posAngle = ConstrainedLinearControlContinuousSystem(A, E, X, B*U)

    # transition posAngle -> deadzone
    X_l3l2 = HPolyhedron([HalfSpace([1.; z], α)])  # x <= 0.03
    r4 = ConstrainedLinearDiscreteSystem(A_trans, X_l3l2)

    # modes
    modes = [m_negAngle, m_deadzone, m_posAngle]

    # reset maps
    resets = [r1, r2, r3, r4]

    # ===========
    # switchings
    # ===========
    switchings = [HybridSystems.AutonomousSwitching()]

    # instantiate hybrid system
    HS = HybridSystem(automaton, modes, resets, switchings)

    return HS
end
