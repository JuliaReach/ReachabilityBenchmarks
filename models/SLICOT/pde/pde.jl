#=
Model: pde.jl
=#
include("../../src/Reachability.jl")

using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen("pde.mat")
    A = read(file, "A")*1.0 # this model provides a matrix with Int components

    # initial set
    n = size(A, 1)
    center0 = zeros(n)
    radius0 = zeros(n)
    center0[65:80] = 0.00125
    center0[81:84] = -0.00175
    radius0[65:84] = 0.00025
    X0 = Hyperrectangle(center0, radius0)

    # input set
    B = read(file, "B")
    U = B * BallInf([0.75], .25)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property(read(matopen("out.mat"), "M")[1,:], 12.), # y < 12
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.003, # time step
#       :blocks => [1],
#       :projection_matrix => sparse(read(matopen("out.mat"), "M")),
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
