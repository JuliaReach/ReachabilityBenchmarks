#=
Model: MNA_5 (10913 variables, 9 inputs)
=#
using Reachability, MAT, SparseArrays

mna5(o::Pair{Symbol, <:Any}...) = mna5(Options(Dict{Symbol, Any}(o)))

function mna5(input_options::Options)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "mna5.mat")
    A = sparse(read(file, "A"))

    # initial set:
    X0 = Hyperrectangle([fill(0.000225, 10); zeros(10903)],
                        [fill(0.000025, 10); zeros(10903)])

    # input set
    B = sparse(19:27, 1:9, fill(-1., 9), size(A, 1), 9)
    U = B * Hyperrectangle([fill(0.1, 5); fill(0.2, 4)], zeros(9))

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # property: x1 < 0.2 && x2 < 0.15
    property = SafeStatesProperty(HPolyhedron([
        HalfSpace(sparsevec([1], [1.0], 10913), 0.2),
        HalfSpace(sparsevec([2], [1.0], 10913), 0.15)]))

    # =======================
    # Problem default options
    # =======================
    # 2D blocks except for the last block, which is 1D
    partition = vcat([(2*i-1:2*i) for i in 1:5456], [10913:10913])

    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [1],
                                  :partition => partition,
                                  :plot_vars => [0, 1],
                                  :assume_sparse => true,
                                  :lazy_expm => true)

    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => [1, 2],
                                  :partition => partition,
                                  :property => property,
                                  :assume_sparse => true,
                                  :lazy_expm => true)
    end

    return (S, merge(problem_options, input_options))
end
