using LazySets

function fiveDimSys_specification()
    # initial set
    X0 = BallInf([1.0, 0.0, 0.0, 0.0, 0.0], 0.1)

    # time horizon: 5 time units
    time_horizon = 5.0

    # specification
    O = Dict{Symbol,Any}(:T => time_horizon)

    return X0, O
end
