#=
Model:  fiveDimSys.jl

This is a five-dimensional model taken from [Girard05]. This also appears as
Example 4.1. page 57 in the thesis [LeGuernic09].

See also Fig. 5.1. page 71 in the same thesis.

[Girard05]       -- Girard, Antoine. "Reachability of uncertain linear systems using zonotopes."
                    International Workshop on Hybrid Systems: Computation and Control.
                    Springer, Berlin, Heidelberg, 2005.
[LeGuernic09]    -- Le Guernic, Colas. Reachability analysis of hybrid systems
                    with linear continuous dynamics. Diss.
                    Université Joseph-Fourier-Grenoble I, 2009.
=#
using Reachability, Plots

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    D = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2.]
    P = [0.6 -0.1 0.1 0.7 -0.2; -0.5 0.7 -0.1 -0.8 0; 0.9 -0.5 0.3 -0.6 0.1;
        0.5 -0.7 0.5 0.6 0.3; 0.8 0.7 0.6 -0.3 0.2]
    A = P * D * inv(P)

    # initial set
    X0 = BallInf([1.0, 0.0, 0.0, 0.0, 0.0], 0.1)

    # input set
    U = Ball2(zeros(5), 0.01)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U) # @system x' = A*x + u, u ∈ U

    # ===============
    # Problem solving
    # ===============
    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :blocks => [1],
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

compute(:N => 100, :T => 1.0); # warm-up
compute(:δ => 0.01, :T => 5.0); # benchmark settings (long)
