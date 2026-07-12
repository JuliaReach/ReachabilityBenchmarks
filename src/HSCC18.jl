using ReachabilityModels, SparseArrays, MAT

M_pde = vec(matread("out_pde.mat")["M"])
M_iss = vec(matread("out_iss.mat")["M"])

y = [false, false, true, false, true, false]
models = ["motor", "building", "pde", "heat", "iss", "beam"]
dirs = [5, 25, M_pde, 133, M_iss, 89]
conditions = [x -> isdisjoint(x, Interval(0.45, 0.6)),
              x -> isdisjoint(x, HalfSpace([-1.], -6e-3)),
              x -> isdisjoint(x, HalfSpace([-1.], -12.)),
              x -> isdisjoint(x, HalfSpace([-1.], -0.1)),
              x -> x ⊆ (Interval(-7., 7.) * 1e-4),
              x -> isdisjoint(x, HalfSpace([-1.], -2100.))]
sols = []
ints = []

for i=1:length(models)
    println(models[i])
    m = fetch_model(models[i])
    meta = fetch_meta(models[i])
    dim = meta["info"]["dim"]
    if !y[i]
        x = [zeros(dirs[i]-1); 1.0; zeros(dim-dirs[i])]
    else
        x = dirs[i]
    end
    @time sol = solve(m, T=20, alg=LGG09(δ=0.001, template=[x, -x]))
    push!(sols, sol)
    if !y[i]
        int = Interval(-ρ(sparsevec([dirs[i]], [-1.0], dim), sols[i]), ρ(sparsevec([dirs[i]], [1.0], dim), sols[i]))
    else
        int = Interval(-ρ(x, sols[i]), ρ(-x, sols[i]))
    end
    push!(ints, int)
    if conditions[i](int)
        println("verified")
    else
        println("not verified")
    end
end
