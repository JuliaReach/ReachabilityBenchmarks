using LazySets

function spikingNeuron_specification()
    # initial set in mode 1
    X0_m1 = Hyperrectangle(; low=[-65.0, -0.2], high=[-60.0, 0.2])
    X0 = [(1, X0_m1)]

    # time horizon: 100 time units
    time_horizon = 100.0

    # specification
    O = Dict{Symbol,Any}(:T => time_horizon)

    return X0, O
end
