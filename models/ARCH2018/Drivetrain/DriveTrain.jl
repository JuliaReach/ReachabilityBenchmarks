function drivetrain(nϴ, opD, t, max_jumps)::AbstractSolution
    system_dimension = 2 * nϴ + 7 + 1
    z = zeros(system_dimension-1)


    # Based on: Avoiding Geometric Intersection Operations in Reachability Analysis of Hybrid Systems
    # Matthias Althoff & Bruce H. Krogh
    #constants
    α = 0.03
    τ_eng = 0.1
    b_l = 5.6
    b_m = 0.
    b_i = 1.
    k_s = 10000.
    k_i = 100000.
    J_l = 140.
    J_m = 0.3
    J_i = 0.01
    γ = 12.
    k_P = 0.5
    k_I = 0.5
    k_D = 0.5
    u = 5

    # common flow
    A = zeros(system_dimension, system_dimension)
    A[1,7], A[1,9] = 1./γ, -1.
    A[2,1], A[2,2], A[2,3], A[2,4], A[2,7], A[2,8] = (-(1/τ_eng)*k_I*γ+k_D*(1/J_m)*k_s),-k_D*(1/J_m - 1),
            (1/τ_eng)*k_I*γ, (1/τ_eng)*k_P*γ, -(1/τ_eng)*k_I + k_D*(1/J_m)*b_m, -(1/τ_eng)*k_I
    #TODO add U for x2' : k_D*γ*u-k_D*(1/(J_m*\gamma))k_s*α
    A[3,4] = 1.
    A[5,6] = 1.
    A[6,5], A[6,6], A[6,2*ϴ+6] = -(1/J_l)*k_ϴ, -(1/J_l)*b_l, (1/J_l)*k_ϴ
    A[7,1], A[7,2], A[7,7] = (1/(J_m*γ))*k_s, 1/J_m, -(1/J_m)*b_m
    #TODO add U for x7' : (1/(J_m*\gamma))*k_s*α
    A[9,1], A[9,8], A[9,9], A[9,10] = (1/J_1)*k_s, -(1/J_1)*k_1,-(1/J_1)*b_1, (1/J_1)
    #TODO add U for x_9' : -(1/J_1)*k_s*α
    J_arr = ones((system_dimension - 1 - 7)/2)
    b_arr = ones((system_dimension - 1 - 7)/2)
    k_arr = ones((system_dimension - 1 - 7)/2)
    i = 10
    while i < system_dimension-1
        el = (system_dimension - 8 - i)/2
        A[i,i+1] = 1.
        A[i+1,5], A[i+1,2*τ+4], A[i+1,2*τ+6], A[i+1,2*τ+7] = (1/J_arr[el])*k_arr[ϴ], -(1/J_arr[el])*k_arr[ϴ-1],
        -(1/J_arr[el])*(k_arr[ϴ-1]+k_arr[ϴ]), -(1/J_arr[el])*b_arr[ϴ]
        i += 2
    end

    # Transition graph (automaton)
    a = LightAutomaton(5);

    add_transition!(a, 1, 2, 1);
    add_transition!(a, 2, 3, 1);
    add_transition!(a, 3, 2, 1);
    add_transition!(a, 3, 4, 1);
    add_transition!(a, 4, 3, 1);

    #negAngleInit
    B = [sparsevec(4:4, [u], system_dimension)]
    X = HPolyhedron([HalfSpace([z; 1.], 0.2)]) # t <= 0.2
    m_negAngleInit = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    #negAngle
    B = [sparsevec(4:4, [u], system_dimension)]
    X = HPolyhedron([HalfSpace([1.; z], -0.03)]) # x <= -0.03
    m_negAngle = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    #Deadzone
    B = [sparsevec(4:4, [u], system_dimension)]
    X = HPolyhedron([HalfSpace([-1.; z], -0.03),  # x >= -0.03
              HalfSpace([1.; z], 0.03)])  # x <= 0.03
    m_deadzone = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    #posAngle
    B = [sparsevec(4:4, [u], system_dimension)]
    X = HPolyhedron([HalfSpace([1; z], 0.03)])  # x <= 0.03 (2.1 for numerical issues)
    m_4 = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    m_posAngle = [m_1, m_2, m_3, m_4];

    # common resets
    A_trans = eye(system_dimension)

    # Transition negAngleInit -> negAngle
    X_l1l2 = HPolyhedron(HalfSpace([z,1.;], 0.2)])  # t <= 0.2
    r1 = ConstrainedLinearDiscreteSystem(A_trans_34, X_l3l4);

    # Transition negAngle -> Deadzone
    X_l2l3 = HPolyhedron([HalfSpace([-1.; z], -0.03)])  # x >= -0.03
    r2 = ConstrainedLinearDiscreteSystem(A_trans, X_l4l2);

    # Transition Deadzone -> negAngle
    X_l3l2 = HPolyhedron([HalfSpace([1.; z], -0.03)])  # x <= -0.03
    r3 = ConstrainedLinearDiscreteSystem(A_trans, X_l4l2);

    # Transition Deadzone -> posAngle
    X_l3l4 = HPolyhedron([HalfSpace([-1.; z], 0.03)])  # x >= 0.03
    r4 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l1);

    # Transition posAngle -> Deadzone
    X_l4l3 = HPolyhedron([HalfSpace([1.; z], 0.03)])  # x >= -0.03
    r5 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l1);

    r = [r1,r2,r3,r4,r5]

    # Switchings
    s = [HybridSystems.AutonomousSwitching()];

    HS = HybridSystem(a, m, r, s);

end
