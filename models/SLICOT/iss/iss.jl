#=
Model: iss.jl

This is a 270-variable model of component 1r (Russian service module) of the
International Space Station (ISS).

The corresponding SpaceEx model and configuration file are iss.xml and iss.cfg.
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen(@relpath "iss.mat")
    A = sparse(read(file, "A"))

    # initial set
    X0 = BallInf(zeros(size(A, 1)), .0001) # -0.0001 <= xi <= 0.0001 for all i

    # input set
    #Uraw = BallInf([0.05], .05) * BallInf([0.9], .1) * BallInf([0.95], .05)
    #Uraw = Hyperrectangle([0.05, 0.9, 0.95], [0.05, 0.1, 0.05])
    Uraw = CartesianProductArray([BallInf([0.05], .05), BallInf([0.9], .1), BallInf([0.95], .05)])
    B = read(file, "B")
    # input U is constant
    U = B * Uraw

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property(read(matopen(@relpath "out.mat"), "M")[1,:], 7e-4), # y < 7e-4
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.01, # time step
#       :blocks => [@block_id(182)],
        :blocks => 68:135, # blocks needed for property
        :assume_sparse => true,
#       :projection_matrix => sparse(read(matopen(@relpath "out.mat"), "M")),
        :plot_vars => [0, 182]
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
