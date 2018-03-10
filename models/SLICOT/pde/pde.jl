#=
Model: pde.jl
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "pde.mat")
    A = read(file, "A")*1.0 # this model provides a matrix with Int components

    # initial set
    n = size(A, 1)
    center0 = zeros(n)
    radius0 = zeros(n)
    center0[65:80] = 0.00125
    center0[81:84] = -0.00175
    radius0[65:84] = 0.00025
    X0 = Hyperrectangle(center0, radius0)

    # input set
    B = read(file, "B")
    U = B * BallInf([0.75], .25)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty(read(matopen(@relpath "out.mat"), "M")[1,:], 12.), # y < 12
#       :vars => [1], # variable for single block analysis
        :vars => 1:42, # variables needed for property
        :partition => [(2*i-1:2*i) for i in 1:42], # 2D blocks
#       :projection_matrix => sparse(read(matopen(@relpath "out.mat"), "M")),
        :plot_vars => [0, 1]
        ), Options(input_options...))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        if options[:mode] == "reach"
            println("Plotting...")
            tic()
            plot(result) # TODO output labels
            @eval(savefig("motor.jl"))
            toc()
        end
    end
end # function

compute(:δ => 0.003, :N => 3); # warm-up
compute(:δ => 0.003, :T => 20.0); # benchmark settings (long)
