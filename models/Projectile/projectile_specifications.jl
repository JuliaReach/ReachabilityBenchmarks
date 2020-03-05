using LazySets

function projectile_specification()
    # initial set
    X0 = Singleton([0.0, 5.0, 100.0, 0.0])

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Dict{Symbol, Any}(:T => time_horizon)

    return X0, O
end
