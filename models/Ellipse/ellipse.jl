#=
Simple 2D ellipse with convex combination initial state.
=#
using Reachability, Plots

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    A = sparse([1, 1, 2, 2], [1, 2, 1, 2], [3., -9., 4., -3.])

    # initial set
    X0_1 = Singleton([1.0, 0.0])
    X0_2 = Singleton([1.5, 0.0])
    alpha = 0.8
    Z0 = Singleton(alpha * X0_1.element + (1.-alpha) * X0_2.element)

    # instantiate continuous LTI system
    S1 = ContinuousSystem(A, X0_1)
    S2 = ContinuousSystem(A, X0_2)
    S3 = ContinuousSystem(A, Z0)

    # ===============
    # Problem solving
    # ===============
    # define solver-specific options
    options = Options(
        :δ => 0.01,
        :approx_model => "nobloating",
        :plot_vars => [1, 2]
        )

    t1 = 1.2
    t2 = 1.0
    t3 = 1.09

    o1 = merge(options, Options(:T => t1))
    o2 = merge(options, Options(:T => t2))
    o3 = merge(options, Options(:T => t3))

    # solve
    result1 = solve(S1, o1)
    result2 = solve(S2, o2)
    result3 = solve(S3, o3)

    # ========
    # Plotting
    # ========
    plot(result1,color="red")
    plot!(result2,color="blue")
    plot!(result3,color="green")
    @eval(savefig(@filename_to_png))
end # function

compute()
