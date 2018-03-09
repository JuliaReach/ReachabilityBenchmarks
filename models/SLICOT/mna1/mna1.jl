#=
Model: MNA_1.jl

This is a 578-variable model.
=#
using Reachability, LazySets, MAT, Plots

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "mna1.mat")
    A = sparse(read(file, "A"))

    # initial set
    X0 = Hyperrectangle([fill(0.00125, 2); zeros(576)], [fill(0.00025, 2); zeros(576)])

    # input set
    B = sparse(570:578, 1:9, fill(-1.0, 9), size(A, 1), 9)
    U = B * Hyperrectangle([fill(0.1, 5); fill(0.2, 4)], zeros(9))

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty([1., 0.], 0.5), # x1 < 0.5
        :vars => [1], # variable needed for property
        :partition=> [(2*i-1:2*i) for i in 1:289], # 2D blocks
        :plot_vars => [0, 1]
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


compute(:δ => 0.001, :N => 3); # warm-up
compute(:δ => 0.001, :T => 20.0); # benchmark settings (long)
