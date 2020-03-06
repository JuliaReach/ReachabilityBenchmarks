using LazySets, MathematicalPredicates

function crane_specification()
    # initial set
    X0 = Hyperrectangle([2.5, 0.0, 0.0, 0.0, 0.0, 0.0],
                        [2.5, 0.0, 0.2, 0.1, 0.0, 0.0])

    # safety property: x1 ≥ -1.8
    property = is_contained_in(HalfSpace([-1., 0.], 1.8))

    # time horizon: 15 time units
    time_horizon = 15.0

    # specification
    O = Dict(:T => time_horizon, :property => property)

    return X0, O
end
