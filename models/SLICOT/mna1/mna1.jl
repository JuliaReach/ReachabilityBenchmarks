#=
Model: MNA_1.jl

This is a 578-variable model.
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen("mna1.mat")
    A = sparse(read(file, "A"))

    # initial set
    X0 = Hyperrectangle([fill(0.00125, 2); zeros(576)], [fill(0.00025, 2); zeros(576)])

    # input set
    B = sparse(570:578, 1:9, fill(-1.0, 9), size(A, 1), 9)
    U = B * Hyperrectangle([fill(0.1, 5); fill(0.2, 4)], zeros(9))

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property([1., 0.], 0.5), # x1 < 0.5
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.001, # time step
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
