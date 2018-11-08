function drivetrain(nϴ, opD, t, max_jumps)::AbstractSolution
    system_dimension = 2 * nϴ + 7 + 1
    z = zeros(system_dimension-1)


    # Based on: Avoiding Geometric Intersection Operations in Reachability Analysis of Hybrid Systems
    # Matthias Althoff & Bruce H. Krogh

    #constants
    constants = Dict(:α => 0.03, :τ_eng => 0.1, :b_l => 5.6, :b_m => 0., :b_i => 1.,
        :k_s => 10000., :k_i => 100000., :J_l => 140., :J_m => 0.3, :J_i => 0.01, :γ => 12.,
            :k_P => 0.5, :k_I => 0.5, :k_D => 0.5, :u => 5)

    # Transition graph (automaton)
    a = LightAutomaton(3);

    add_transition!(a, 1, 2, 1);
    add_transition!(a, 2, 1, 1);
    add_transition!(a, 2, 3, 1);
    add_transition!(a, 3, 2, 1);

    # common U
    U = Singleton([1.0]);

    #negAngle
    constants[:k_s] = 10000.
    constants[:α] = -0.03
    A = get_dynamics(constants, system_dimension)
    B = get_b(constants)
    X = HPolyhedron([HalfSpace([1.; z], constants[:α])]) # x <= -α
    m_negAngle = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);


    #Deadzone
    constants[:k_s] = 0
    constants[:α] = -0.03
    A = get_dynamics(constants, system_dimension)
    B = get_b(constants)
    X = HPolyhedron([HalfSpace([-1.; z], constants[:α]),  # x >= -⁠α
              HalfSpace([1.; z], -constants[:α])])  # x <= 0.03
    m_deadzone = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    #posAngle
    constants[:k_s] = 10000.
    constants[:α] = 0.03
    A = get_dynamics(constants, system_dimension)
    B = get_b(constants)
    X = HPolyhedron([HalfSpace([1; z], constants[:α])])  # x <= α (2.1 for numerical issues)
    m_4 = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

    m_posAngle = [m_1, m_2, m_3, m_4];

    # common resets
    A_trans = eye(system_dimension)

    # Transition negAngle -> Deadzone
    X_l1l2 = HPolyhedron([HalfSpace([-1.; z], -constants[:α])])  # x >= -0.03
    r2 = ConstrainedLinearDiscreteSystem(A_trans, X_l4l2);

    # Transition Deadzone -> negAngle
    X_l2l1 = HPolyhedron([HalfSpace([1.; z], -constants[:α])])  # x <= -0.03
    r3 = ConstrainedLinearDiscreteSystem(A_trans, X_l4l2);

    # Transition Deadzone -> posAngle
    X_l2l3 = HPolyhedron([HalfSpace([-1.; z], constants[:α])])  # x >= 0.03
    r4 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l1);

    # Transition posAngle -> Deadzone
    X_l3l2 = HPolyhedron([HalfSpace([1.; z], constants[:α])])  # x >= -0.03
    r5 = ConstrainedLinearDiscreteSystem(A_trans, X_l2l1);

    r = [r1,r2,r3,r4,r5]

    # Switchings
    s = [HybridSystems.AutonomousSwitching()];

    HS = HybridSystem(a, m, r, s);

    
end



function get_dynamics(constants::Dictionary, system_dimension::Int)
    α = constants[:α]
    τ_eng = constants[:τ_eng]
    b_l = constants[:b_l]
    b_m = constants[:b_m]
    b_i = constants[:b_i]
    k_s = constants[:k_s]
    k_i = constants[:k_i]
    J_l = constants[:J_l]
    J_m = constants[:J_m]
    J_i = constants[:J_i]
    γ = constants[:γ]
    k_P = constants[:k_P]
    k_I = constants[:k_I]
    k_D = constants[:k_D]
    u = constants[:u]

    # common flow
    A = zeros(system_dimension, system_dimension)

    A[1,7] = 1.0/γ
    A[1,9] =  -1.

    A[2,1] = -(1.0/τ_eng)*k_I*γ+k_D*(1.0/(J_m*γ*τ_eng)*k_s
    A[2,2] = -k_D*(1.0/J_m - 1)/τ_eng
    A[2,3] = (1.0/τ_eng)*k_I*γ
    A[2,4] = (1.0/τ_eng)*k_P*γ
    A[2,7] = (-k_P + k_D*(1.0/J_m)*b_m)/τ_eng
    A[2,8] = -(1.0/τ_eng)*k_I

    A[3,4] = 1.

    A[5,6] = 1.

    A[6,5] = -(1.0/J_l)*k_ϴ
    A[6,6] =  -(1.0/J_l)*b_l
    A[6,2 * ϴ+6] = (1.0/J_l)*k_ϴ

    A[7,1] = -(1.0/(J_m*γ))*k_s
    A[7,2] = 1.0/J_m
    A[7,7] = -(1.0/J_m)*b_m

    J_arr = fill(J_i, Int((system_dimension - 1 - 7)/2))
    b_arr = fill(b_i, Int((system_dimension - 1 - 7)/2))
    k_arr = fill(k_i, Int((system_dimension - 1 - 7)/2))

    A[9,1] = (1.0/J_arr[1])*k_s
    A[9,8] = -(1.0/J_arr[1])*k_arr[1]
    A[9,9] = -(1.0/J_arr[1])*b_arr[1]
    A[9,10] = (1.0/J_arr[1])*k_arr[1]

    i = 10
    while i < system_dimension-1
        el = (system_dimension - 8 - i)/2
        A[i,i+1] = 1.
        A[i+1,5] = (1.0/J_arr[el])*k_arr[el]
        A[i+1,2*el+4] = (1.0/J_arr[el])*k_arr[el-1]
        A[i+1,2*el+6] = -(1.0/J_arr[el])*(k_arr[el-1]+k_arr[el])
        A[i+1,2*el+7] =  -(1.0/J_arr[el])*b_arr[el]
        i += 2
    end
end

function get_b(constants::Dictionary)
    x2_u = (constants[:k_D]*constants[:γ]*constants[:u]-constants[:k_D]*
            (1.0/(constants[:J_m]*constants[:γ]))*constants[:k_s]*constants[:α])/constants[:τ_eng]
    x4_u = constants[:u]
    x7_u = (1.0/(constants[:J_m]*constants[:γ]))*constants[:k_s]*constants[:α]
    x9_u = -(1.0/constants[:J_i])*constants[:k_s]*constants[:α]

    return sparsevec([2, 4, 7, 9], [x2_u, x4_u, x7_u, x9_u], system_dimension)
end
