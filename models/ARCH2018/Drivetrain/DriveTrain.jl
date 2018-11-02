function drivetrain(nϴ, opD, t, max_jumps)::AbstractSolution
    system_dimension = 2 * nϴ + 7 + 1
    z = zeros(system_dimension-1)

    # common flow
    A = zeros(system_dimension, system_dimension)
    A[1,7], A[1,9] = 1./12., -1.
    A[3,4] = 1.
    A[5,6] = 1.
    A[system_dimension,system_dimension] = 1

    i = 8
    while i < system_dimension-1
        A[i,i+1] = 1.
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
    B = [sparsevec(4:4, [-5], system_dimension)]
    X = HPolyhedron([HalfSpace([1; z], 0.2)]) # t <= 0.2
    m_negAngleInit = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    #negAngle
    B = [sparsevec(4:4, [5], system_dimension)]
    X = HPolyhedron([HalfSpace([1; z], -0.03)]) # x <= -0.03
    m_negAngle = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    #Deadzone
    B = [sparsevec(4:4, [5], system_dimension)]
    X = HPolyhedron([HalfSpace([-1; z], -0.03),  # x >= -0.03
              HalfSpace([1; z], 0.03)])  # x <= 0.03
    m_deadzone = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    #posAngle
    B = [sparsevec(4:4, [5], system_dimension)]
    X = HPolyhedron([HalfSpace([1; z], 0.03)])  # x <= 0.03 (2.1 for numerical issues)
    m_4 = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    m_posAngle = [m_1, m_2, m_3, m_4];

    # common resets
    A_trans = eye(system_dimension)

    # Transition negAngleInit -> negAngle
    X_l1l2 = HPolyhedron(HalfSpace([z,1.;], 0.2)])  # t <= 0.2
    A_trans_34 = eye(system_dimension)
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
