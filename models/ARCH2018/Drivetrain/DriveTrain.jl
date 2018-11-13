using HybridSystems, LazySets
if VERSION >= v"0.7"
    using LinearAlgebra, SparseArrays
end

function drivetrain(nϴ)::HybridSystem
    system_dimension = 2 * nϴ + 7
    z = zeros(system_dimension)

    # Based on: M. Althoff & B. H. Krogh - Avoiding Geometric Intersection
    # Operations in Reachability Analysis of Hybrid Systems

    # constants
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

    function get_dynamics(k_s, α)
        # common flow
        A = spzeros(system_dimension, system_dimension)

        J_arr = fill(J_i, nϴ)
        b_arr = fill(b_i, nϴ)
        k_arr = fill(k_i, nϴ)

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

        A[6,5] = -(1.0/J_l) * k_arr[nϴ]
        A[6,6] =  -(1.0/J_l) * b_l
        A[6,2 * nϴ+6] = (1.0/J_l) * k_arr[nϴ]

        A[7,1] = -(1.0/(J_m*γ))*k_s
        A[7,2] = 1.0/J_m
        A[7,7] = -(1.0/J_m)*b_m

        i = 10
        if (nϴ > 1)
            A[8, 9] = 1.
            A[9,1] = (1.0/J_arr[1])*k_s
            A[9,8] = -(1.0/J_arr[1])*k_arr[1]
            A[9,9] = -(1.0/J_arr[1])*b_arr[1]
            A[9,10] = (1.0/J_arr[1])*k_arr[1] # ??? maybe it is a typo in the paper?
        else
            i = 8
        end
        while i < system_dimension
            el = (system_dimension - i)/2
            A[i,i+1] = 1.
            A[i+1,5] = (1.0/J_arr[el])*k_arr[el]
            A[i+1,2*el+4] = (1.0/J_arr[el])*k_arr[el-1]
            A[i+1,2*el+6] = -(1.0/J_arr[el])*(k_arr[el-1]+k_arr[el])
            A[i+1,2*el+7] =  -(1.0/J_arr[el])*b_arr[el]
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

        return sparsevec([2, 4, 7, 9], [x2_u, x4_u, x7_u, x9_u], system_dimension)
    end


    # transition graph (automaton)
    a = LightAutomaton(3);

    add_transition!(a, 1, 2, 1);
    add_transition!(a, 2, 1, 1);
    add_transition!(a, 2, 3, 1);
    add_transition!(a, 3, 2, 1);

    # common U
    U = Singleton([1.0]);
    # common resets
    A_trans = eye(system_dimension)

    # negAngle
    A = get_dynamics(k_s, α_neg)
    B = get_b(k_s, α_neg)
    X = HPolyhedron([HalfSpace([1.; z], α_neg)]) # x <= -α
    m_negAngle = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    # transition negAngle -> deadzone
    X_l1l2 = HPolyhedron([HalfSpace([-1.; z], α)])  # x >= -0.03
    r1 = ConstrainedLinearDiscreteSystem(A_trans, X_l1l2);


    # deadzone
    A = get_dynamics(k_s_zero, α_neg)
    B = get_b(k_s_zero, α_neg)
    X = HPolyhedron([HalfSpace([-1.; z], α),  # x >= -α
              HalfSpace([1.; z], α)])  # x <= 0.03
    m_deadzone = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);
    # transition deadzone -> negAngle
    X_l2l1 = HPolyhedron([HalfSpace([1.; z], α_neg)])  # x <= -0.03
    r2 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l1);
    # transition deadzone -> posAngle
    X_l2l3 = HPolyhedron([HalfSpace([-1.; z], α_neg)])  # x >= 0.03
    r3 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l3);


    #posAngle
    A = get_dynamics(k_s, α)
    B = get_b(k_s, α)
    X = HPolyhedron([HalfSpace([-1.; z], α_neg)])  # x >= α (2.1 for numerical issues)
    m_posAngle = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);
    # transition posAngle -> deadzone
    X_l3l2 = HPolyhedron([HalfSpace([1.; z], α)])  # x <= 0.03
    r4 = ConstrainedLinearDiscreteSystem(A_trans, X_l3l2);

    modes = [m_negAngle, m_deadzone, m_posAngle];

    r = [r1,r2,r3,r4];

    # switchings
    s = [HybridSystems.AutonomousSwitching()];

    HS = HybridSystem(a, m, r, s);

    return HS;
end
