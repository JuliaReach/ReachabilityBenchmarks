#=
Model: Beam (348 variables, 1 input)
=#
using Reachability, MAT, SparseArrays


#get_filename() = @relpath "beam.mat"

function beam(input_options=nothing; filename=@relpath "beam.mat")

    # =====================
    # Problem specification
    # =====================
    file = matopen(filename)
    A = read(file, "A")

    # initial set
    # - x1-x300 are 0,
    # - the rest is in [0.002, 0.0015]
    X0 = Hyperrectangle([zeros(300); fill(0.00175, 48)],
                        [zeros(300); fill(0.00025, 48)])

    # input set
    B = read(file, "B")
    U = BallInf([0.5], 0.3)

    # instantiate continuous LTI system
    P = IVP(CLCCS(A, B, nothing, U), X0)

    # safety property: x89 < 2100
    property = SafeStatesProperty(HalfSpace(sparsevec([89], [1.0], 348), 2100.))

    # =======================
    # Problem options
    # =======================
    O = Options(:property => property, :plot_vars => [0, 89], :T=>20.0)

    # algorithm-specific options
    partition = [(2*i-1:2*i) for i in 1:174] # 2D blocks
    O_BFFPSV18 = Options(:vars => [89], :partition => partition, :δ=>1e-3)
    if input_options != nothing
        O_BFFPSV18 = merge(problem_options, input_options)
    end
    return P, O, BFFPSV18(O_BFFPSV18)
end
