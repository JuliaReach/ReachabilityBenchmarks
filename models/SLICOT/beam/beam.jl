#=
Model: beam.jl

This is a 348-variable model.
=#
using Reachability, MAT, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)

    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "beam.mat")
    A = sparse(read(file, "A"))

    # initial set
    # - x1-x300 are 0.0,
    # - the rest is in [0.002, 0.0015]
    X0 = Hyperrectangle([zeros(300); fill(0.00175, 48)], [zeros(300); fill(0.00025, 48)])

    # input set
    B = sparse(read(file, "B"))
    U = B * BallInf([0.5], 0.3)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [89],
                                 :partition => [(2*i-1:2*i) for i in 1:174],
                                 :plot_vars => [0, 89])

    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => [89],
                                  :partition => [(2*i-1:2*i) for i in 1:174],
                                  :property => LinearConstraintProperty(sparsevec([89], [1.0], 348), 2100.)) # x89 < 2100
    end

    result = solve(S, merge(input_options, problem_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@relpath "beam.png"))
        toc()
    end
end # function

# Reach tube computation in dense time
compute(:δ => 1e-3, :N => 3, :mode=>"reach", :verbosity => "info"); # warm-up
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach", :verbosity => "info"); # benchmark settings (long)
