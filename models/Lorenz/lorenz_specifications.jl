using LazySets

function lorenz_specification()
    # initial set
    X0 = Hyperrectangle(; low=[0.9, 0.0, 0.0], high=[1.1, 0.0, 0.0])

    # time horizon: 10 time units
    time_horizon = 10.0

    # specification
    O = Dict{Symbol,Any}(:T => time_horizon)

    return X0, O
end
