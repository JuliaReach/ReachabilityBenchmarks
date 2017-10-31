#=
Model : crane.jl

=#
using Reachability, LazySets

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    A = sparse([0. 1. 0. 0. 0. 0. ;
        -0.417533 -3.1931759963 39.24 0. -14.825331 11.123344 ;
        0. 0. 0. 1. 0. 0. ;
        0.0417533 0.31931759963 -4.905 0. 1.4825331 -1.1123344 ;
        0.0638407957 -0.32473339016573 0. 0. -3.7332068901 -0.7007592976 ;
        0.0853437452 -0.72366802635628 0. 0. -5.9714023436 -2.2736115136])

    # initial set
    X0 = CartesianProductArray([BallInf([2.5], 2.5),
        Singleton([0.]), BallInf([0.], 0.2), BallInf([0.], 0.1), Singleton([0.,0.])])

    # input set
    U = VoidSet(6)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property([-1., 0.], 1.8), # -x1 < 1.8 == x1 > -1.8
        :T => 15., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.01, # time step
        :blocks => [1],
        :plot_vars => [0, 1]
        ), Options(input_options...))

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
