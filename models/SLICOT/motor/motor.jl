#=
Model: motor.jl

For the SpaceEx model, see [1] Motor.xml and the configuration file Motor.cfg.
Output variables are x1 and x2.

[1] BakDuggirala2017cavrepeatability/arch_benchmarks/Motor/Motor.xml
=#
include("../../src/Reachability.jl")

using Reachability, LazySets

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    I = [1, 2, 2, 3, 3, 3, 3, 4, 5, 6, 6, 7, 7, 7, 7, 8]
    J = [2, 3, 2, 1, 2, 3, 4, 1, 6, 7, 6, 5, 6, 7, 8, 5]
    vals = [1, 8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0, 1.0, 8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0]
    A = sparse(I, J, vals)

    # initial set
    X0 = Hyperrectangle([0.00225, 0.0, 0.0, 0.0, 0.00125, 0.0, 0.0, 0.0], [0.00025, 0.0, 0.0, 0.0, 0.00025, 0.0, 0.0, 0.0])

    # input set
    B = sparse([4, 8], [1, 2], [-1.0, -1.0])
    U = B * Hyperrectangle([0.23, 0.3], [0.07, 0.1])

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property(Clause([LinearConstraint([1.; zeros(7)], 0.35), LinearConstraint([zeros(4); 1.; zeros(3)], 0.45)])), # x1 < 0.35 || x5 < 0.45
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.001, # time step
        :blocks => [1, 3], # blocks needed for property
        :plot_vars => [0, 5]
        ), Options(Dict{Symbol,Any}(input_options)))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        options_plot = Options(
            :plot_vars => options[:plot_vars],
            :plot_name => @filename_to_png
#           :plot_indices => range_last_x_percent(length(result), 10, 3)
            )
        plot(result, options_plot)
        toc()
    end
end # function
nothing
