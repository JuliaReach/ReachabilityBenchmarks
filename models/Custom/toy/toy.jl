#=
Model: toy.jl

Toy model to test the implementation on certain cases.
=#
using Reachability

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    A = sparse([1, 2], [1, 2], [0.1, 0.1])

    # initial set
    X0 = BallInf([1.0, 1.0], 0.5)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property([1., 0.], 20.), # x1 < 20
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.01, # time step
        :blocks => [1],
        :plot_vars => [0, 1]
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
