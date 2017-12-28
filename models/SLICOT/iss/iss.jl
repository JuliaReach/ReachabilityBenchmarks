#=
Model: iss.jl

This is a 270-variable model of component 1r (Russian service module) of the
International Space Station (ISS).

The corresponding SpaceEx model and configuration file are iss.xml and iss.cfg.
=#
using Reachability, MAT, Plots

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
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

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty(read(matopen(@relpath "out.mat"), "M")[1,:], 7e-4), # y < 7e-4
#       :blocks => [@block_id(182)],
        :blocks => 68:135, # blocks needed for property
        :assume_sparse => true,
#       :projection_matrix => sparse(read(matopen(@relpath "out.mat"), "M")),
        :plot_vars => [0, 182]
        ), Options(input_options...))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        #project_output = options[:projection_matrix] != nothing
        #:plot_labels => add_plot_labels(options[:plot_vars], project_output)
        plot(result) #TODO: labels, project output
        @eval(savefig(@filename_to_png))
        toc()
    end
end # function

compute(:N => 10, :T => 20.0); # warm-up
compute(:δ => 0.01, :T => 20.0); # benchmark settings (long)
