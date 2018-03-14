#=
Model: MNA_1.jl

This is a 578-variable model.
=#
using Reachability, MAT, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)
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
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [1],
                                  :partition => [(2*i-1:2*i) for i in 1:289],
                                  :plot_vars => [0, 1],
                                  :assume_sparse => false)
    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => [1], # variables needed for property
                                  :partition => [(2*i-1:2*i) for i in 1:289],
                                  :property => LinearConstraintProperty(sparsevec([1], [1.0], 578), 0.5), # x1 < 0.5
                                  :assume_sparse => false)
    end

    result = solve(S, merge(input_options, problem_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@relpath "mna1.png"))
        toc()
    end
end # function


# Reach tube computation in dense time
compute(:δ => 1e-3, :N => 3, :mode=>"reach", :verbosity => "info"); # warm-up
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach", :verbosity => "info"); # benchmark settings (long)
