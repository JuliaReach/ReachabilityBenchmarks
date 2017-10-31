#=
Model: heat.jl

This is a 200-variable model.
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

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

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property([1., 0.], 0.1), # x133 < 0.1
        :blocks => [@block_id(133)],
        :plot_vars => [0, 133]
        ), Options(input_options...))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result) # TODO: project_output
        @eval(savefig(@filename_to_png))
        toc()
    end
end # function

compute(:N => 10, :T => 20.0); # warm-up
compute(:δ => 0.001, :T => 20.0); # benchmark settings (long)
