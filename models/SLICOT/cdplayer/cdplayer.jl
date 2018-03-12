#=
Model: cdplayer.jl

This is a 120-variable model.
=#
using Reachability, MAT, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "cdplayer.mat")
    A = read(file, "A")

    # initial set
    X0 = BallInf(zeros(120), 1.0)

    # input set
    B = read(file, "B")
    U = B * BallInf([0.0, 0.0], 1.0)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [1],
                                  :partition => [(2*i-1:2*i) for i in 1:60],
                                  :plot_vars => [0, 1],
                                  :assume_sparse => true)
    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => [1, 2],
                                  :partition => [(2*i-1:2*i) for i in 1:60],
                                  :property => LinearConstraintProperty(sparsevec([1, 2], [2., -3.], 120), 450.8), # 2*x1 -3*x2 < 450.8
                                  :assume_sparse => true)
    end

    result = solve(S, merge(input_options, problem_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@relpath "cdplayer.png"))
        toc()
    end
end # function

# Reach tube computation in dense time.
# warning: this model behaves badly with :δ => 1e-3, try smaller δ.
compute(:δ => 1e-3, :N => 3, :mode=>"reach", :verbosity => "info"); # warm-up
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach", :verbosity => "info"); # benchmark settings (long)
