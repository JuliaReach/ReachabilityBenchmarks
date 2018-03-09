#=
Model: cdplayer.jl

This is a 120-variable model.
=#
using Reachability, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "cdplayer.mat")
    A = sparse(read(file, "A"))

    # initial set
    X0 = BallInf(zeros(120), 1.0)

    # input set
    B = sparse(read(file, "B"))
    U = B * BallInf([0.0], 1.0)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty([2., -3.], 450.8), # 2*x1 -3*x2 < 450.8
#       :vars => [1], # variable for single block analysis
        :vars => [1, 2], # variables needed for property
        :partition=> [(2*i-1:2*i) for i in 1:60], # 2D blocks
        :assume_sparse => true,
        :plot_vars => [0, 1]
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

compute(:δ => 0.0001, :N => 3); # warm-up
compute(:δ => 0.0001, :T => 20.0); # benchmark settings (long)

