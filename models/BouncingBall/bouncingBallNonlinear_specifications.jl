using LazySets

function bouncingBallNonlinear_specification()
    # initial set in the "down" mode
    X0_down = Hyperrectangle(; low=[4.9, -0.2], high=[5.1, 0.0])
    X0 = [(2, X0_down)]

    # time horizon: 6 time units
    time_horizon = 6.0

    # specification
    O = Dict{Symbol,Any}(:T => time_horizon)

    return X0, O
end
