using LazySets, Reachability

function motor_specification()
    # initial set:
    # xᵢ ∈ [0.002, 0.0025] for i = 1
    # xᵢ ∈ [0.001, 0.0015] for i = 5
    # xᵢ = 0 otherwise
    X0 = Hyperrectangle(low=[0.002, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0],
                        high=[0.0025, 0.0, 0.0, 0.0, 0.0015, 0.0, 0.0, 0.0])
    X1 = Hyperrectangle([0.00225, 0.0, 0.0, 0.0, 0.00125, 0.0, 0.0, 0.0],
                        [0.00025, 0.0, 0.0, 0.0, 0.00025, 0.0, 0.0, 0.0])

    # safety property: x1 ≤ 0.35 ∨ x5 ≤ 0.45
    property = Disjunction([
        SafeStatesProperty(HalfSpace([1.; zeros(7)], 0.35)),
        SafeStatesProperty(HalfSpace([zeros(4); 1.; zeros(3)], 0.45))])

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Options(:T => time_horizon, :property => property)

    return X0, O
end
