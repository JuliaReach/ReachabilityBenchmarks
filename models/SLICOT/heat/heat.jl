#=
Model: heat.jl

This is a 200-variable model.
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen("heat.mat")
    A = sparse(read(file, "A"))

    # initial set
    # - x1-x300 are 0.0,
    # - the rest is in [0.002, 0.0015]
    X0 = Hyperrectangle([fill(0.6125, 2); zeros(198)], [fill(0.0125, 2); zeros(198)])

    # input set
    B = sparse([67], [1], [1.0], size(A, 1), 1)
    U = B * BallInf([0.0], 0.5)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property([1., 0.], 0.1), # x133 < 0.1
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.001, # time step
        :blocks => [@block_id(133)],
        :plot_vars => [0, 133]
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
