#=
Model: toy.jl

Toy model to test the implementation on certain cases.
=#
using Reachability

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    A = sparse([1, 2], [1, 2], [0.1, 0.1])

    # initial set
    X0 = BallInf([1.0, 1.0], 0.5)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0)

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property([1., 0.], 20.), # x1 < 20
        :blocks => [1],
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
    en
end # function

compute(:N => 10, :T => 20.0); # warm-up
compute(:δ => 0.01, :T => 20.0); # benchmark settings (long)
