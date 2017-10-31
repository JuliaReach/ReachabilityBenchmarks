#=
Model: fom.jl

This is a 1006 x 1006 dimensional model with 1 input.
=#
using Reachability, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen(@relpath "fom.mat")
    A = sparse(read(file, "A"))*1.0  # this model provides a matrix with Int components

    # initial set
    X0 = Hyperrectangle(zeros(1006), [fill(0.0001, 400); zeros(606)])

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
        :property => Property(read(matopen(@relpath "out.mat"), "M")[1,:], 185.), # y < 185
#       :blocks => [1],
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
        plot(result) # TODO: project_output
        @eval(savefig(@filename_to_png))
        toc()
    end
end # function

compute(:N => 10, :T => 20.0); # warm-up
compute(:δ => 0.05, :T => 20.0); # benchmark settings (long)
