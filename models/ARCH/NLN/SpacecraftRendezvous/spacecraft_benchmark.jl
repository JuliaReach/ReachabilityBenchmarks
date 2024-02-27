𝑃, 𝑂 = spacecraft_rendezvous(; T=200.0, orderT=10, orderQ=2, abs_tol=1e-28, max_steps=5000);

𝑂jets = Options(:orderT => orderT, :orderQ => orderQ, :abs_tol => abs_tol, :max_steps => max_steps)

sol = solve(𝑃, 𝑂, TMJets(𝑂jets), LazyDiscretePost(:check_invariant_intersection => true))

# first mode
#p = IVP(𝑃.s.modes[1], 𝑃.x0[1][2])
#solve(p, 𝑂, op=TMJets(𝑂jets))
# compute projection onto the plot variables (project_reachset doesn't currently
# work for a limitation in the hybrid solve, which does not preserve the type of opC)
