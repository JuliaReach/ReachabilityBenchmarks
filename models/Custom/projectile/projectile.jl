#=
Model: projectile.jl
=#
include("../../src/Reachability.jl")

using Reachability, LazySets

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    A = sparse([0. 0.5 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.7 ; 0. 0. 0. 0.])

    # initial set
    X0 = Singleton([0.,5.,100.,0])

    # input set
    U = Singleton([0.,0.,0.,-9.81])

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.5, # time step
        :blocks => [1, 2],
        :plot_vars => [1, 3]
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
