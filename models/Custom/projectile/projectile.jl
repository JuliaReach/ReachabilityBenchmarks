#=
Model: projectile.jl
=#
using Reachability, MathematicalSystems

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    # dynamics
    A = [0. 0.5 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.7 ; 0. 0. 0. 0.]
    b = [0., 0., 0., -9.81]

    # initial set
    X0 = Singleton([0.,5.,100.,0])

    # instantiate continuous LTI system
    S = IVP(AffineContinuousSystem(A, b), X0)

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
    end
end # function

compute(:N => 10, :T => 20.0); # warm-up
compute(:δ => 0.5, :T => 20.0); # benchmark settings (long)
