#=
Model: fom.jl

This is a 1006 x 1006 dimensional model with 1 input.
=#
using Reachability, MAT, Plots

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "fom.mat")
    A = sparse(read(file, "A"))*1.0  # this model provides a matrix with Int components

    # initial set
    X0 = Hyperrectangle(zeros(1006), [fill(0.0001, 400); zeros(606)])

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
        :property => LinearConstraintProperty(read(matopen(@relpath "out.mat"), "M")[1,:], 185.), # y < 185
#       :vars => [1], # variable for single block analysis
        :vars => 1:1006, # variables needed for property
        :partition => [(2*i-1:2*i) for i in 1:503], # 2D blocks
        :assume_sparse => true,
#       :projection_matrix => sparse(read(matopen(@relpath "out.mat"), "M")),
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

compute(:δ => 0.05, :N => 3); # warm-up
compute(:δ => 0.05, :T => 20.0); # benchmark settings (long)
