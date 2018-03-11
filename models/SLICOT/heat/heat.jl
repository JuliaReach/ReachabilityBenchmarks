#=
Model: heat.jl

This is a 200-variable model.
=#
using Reachability, MAT, Plots

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "heat.mat")
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

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty([1., 0.], 0.1), # x133 < 0.1
        :vars => [133], # variable needed for property
        :partition => [(2*i-1:2*i) for i in 1:100], # 2D blocks
        :plot_vars => [0, 133]
        ), Options(input_options...))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@filename_to_png))
        toc()
    end
end # function

compute(:δ => 0.001, :N => 3); # warm-up
compute(:δ => 0.001, :T => 20.0); # benchmark settings (long)
