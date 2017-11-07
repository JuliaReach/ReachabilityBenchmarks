#=
Model: toy.jl

Toy model to test the implementation on certain cases.
=#
using Reachability

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    n = 10
    A = randn(n, n)

    # initial set
    X0 = BallInf(ones(n), 0.1)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0)

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
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
