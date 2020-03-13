using LazySets, MathematicalPredicates

function linearSwitching_specification()
    # initial set in mode 1
    X0_m1 = Singleton([3.1, 4.0, 0.0, 0.0, 0.0])
    X0 = [(1, X0_m1)]

    # safety property: x1 > -1.2
    border = HalfSpace([1.0, 0.0, 0.0, 0.0, 0.0], -1.2)  # x1 ≤ -1.2
    property = is_disjoint_from(border)

    # time horizon: 1 time unit
    time_horizon = 1.0

    # specification
    O = Dict(:T => time_horizon, :property => property)

    return X0, O
end
