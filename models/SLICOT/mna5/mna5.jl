#=
Model: MNA_5.jl

This is a 10913-variable model.
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "mna5.mat")
    A = sparse(read(file, "A"))

    # initial set:
    X0 = Hyperrectangle([fill(0.000225, 10); zeros(10903)], [fill(0.000025, 10); zeros(10903)])

    # input set
    B = sparse(19:27, 1:9, fill(-1., 9), size(A, 1), 9)
    U = B * Hyperrectangle([fill(0.1, 5); fill(0.2, 4)], zeros(9))

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty([Clause([LinearConstraint([1., 0.], 0.2)]), Clause([LinearConstraint([0., 1.], 0.15)])]), # x1 < 0.2 && x2 < 0.15
        :blocks => [1],
        :assume_sparse => true,
        :lazy_expm => true,
        :plot_vars => [0, 1]
        ), Options(input_options...))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        options_plot = Options(
            :plot_vars => options[:plot_vars],
            :plot_name => @filename_to_png
#           :plot_indices => range_last_x_percent(length(result), 10, 3)
            )
        plot(result, options_plot)
        toc()
    end
end # function

compute(:N => 10, :T => 20.0); # warm-up
compute(:δ => 0.05, :T => 20.0); # benchmark settings (long)
