#=
Model: projectile.jl
=#
using Reachability

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    A = [0. 0.5 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.7 ; 0. 0. 0. 0.]

    # initial set
    X0 = Singleton([0.,5.,100.,0])

    # input set
    U = Singleton([0.,0.,0.,-9.81])

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :blocks => [1, 2],
        :plot_vars => [1, 3]
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
compute(:δ => 0.5, :T => 20.0); # benchmark settings (long)
