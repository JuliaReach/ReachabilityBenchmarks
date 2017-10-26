#=
Model: fom.jl

This is a 1006 x 1006 dimensional model with 1 input.
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen(@relpath "fom.mat")
    A = sparse(read(file, "A"))*1.0  # this model provides a matrix with Int components

    # initial set
    X0 = Hyperrectangle(zeros(1006), [fill(0.0001, 400); zeros(606)])

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
        :property => Property(read(matopen(@relpath "out.mat"), "M")[1,:], 185.), # y < 185
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.005, # time step
#       :blocks => [1],
        :assume_sparse => true,
#       :projection_matrix => sparse(read(matopen(@relpath "out.mat"), "M")),
        :plot_vars => [0, 1]
        ), Options(Dict{Symbol,Any}(input_options)))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        project_output = options[:projection_matrix] != nothing
        options_plot = Options(
            :plot_vars => options[:plot_vars],
            :plot_labels => add_plot_labels(options[:plot_vars], project_output),
            :plot_name => @filename_to_png
#           :plot_indices => range_last_x_percent(length(result), 10, 3)
            )
        plot(result, options_plot)
        toc()
    end
end # function
nothing
