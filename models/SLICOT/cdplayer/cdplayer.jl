#=
Model: cdplayer.jl

This is a 120-variable model.
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen("cdplayer.mat")
    A = sparse(read(file, "A"))

    # initial set
    X0 = BallInf(zeros(120), 1.0)

    # input set
    B = sparse(read(file, "B"))
    U = B * BallInf([0.0], 1.0)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property([2., -3.], 450.8), # 2*x1 -3*x2 < 450.8
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.0001, # time step
        :blocks => [1],
        :assume_sparse => true,
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
