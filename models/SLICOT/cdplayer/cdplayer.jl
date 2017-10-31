#=
Model: cdplayer.jl

This is a 120-variable model.
=#
using Reachability, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen(@relpath "cdplayer.mat")
    A = sparse(read(file, "A"))

    # initial set
    X0 = BallInf(zeros(120), 1.0)

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
        :property => Property([2., -3.], 450.8), # 2*x1 -3*x2 < 450.8
        :blocks => [1],
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

compute(:N => 10, :T => 20.0); # warm-up
compute(:δ => 0.0001, :T => 20.0); # benchmark settings (long)

