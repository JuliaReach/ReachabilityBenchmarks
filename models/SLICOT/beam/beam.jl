#=
Model: beam.jl

This is a 348-variable model.
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "beam.mat")
    A = sparse(read(file, "A"))

    # initial set
    # - x1-x300 are 0.0,
    # - the rest is in [0.002, 0.0015]
    X0 = Hyperrectangle([zeros(300); fill(0.00175, 48)], [zeros(300); fill(0.00025, 48)])

    # input set
    B = sparse(read(file, "B"))
    U = B * BallInf([0.5], 0.3)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty([1., 0.], 2100.), # x89 < 2100
        :vars => [89], # variable needed for property
        :partition => [(2*i-1:2*i) for i in 1:174], # 2D blocks
        :plot_vars => [0, 89]
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

compute(:δ => 0.0005, :N => 3); # warm-up
compute(:δ => 0.0005, :T => 20.0); # benchmark settings (long)
