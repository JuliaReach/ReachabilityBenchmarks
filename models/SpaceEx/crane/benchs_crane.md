results crane.jl

(numeros sin contar plot, plot_polygon line is commented)

δ = 0.01
ɛ = 0.01

# perform intermediate overapproximations?
intermediate_overapproximations = true

we do:

N=400
X0 = CartesianProductArray([BallInf([2.5], 2.5), Singleton([0.]), BallInf([0.], 0.2), BallInf([0.], 0.1), Singleton([0.,0.])])

=======================================================
=======================================================
algorithm="no_wrapping_effect_homogeneous"
# instantiate the continuous dynamical system, no input
D = ContinuousSystem(A, X0)

julia> @time compute(400)
System Construction...
elapsed time: 1.8297e-5 seconds
Time Discretization...
elapsed time: 0.000453649 seconds
Reachable LazySets Computation...
no_wrapping_effect_homogeneous
elapsed time: 0.123645074 seconds
Projections...
elapsed time: 0.026079148 seconds
  0.150575 seconds (1.39 M allocations: 64.673 MiB, 23.22% gc time)
0.026079148

BenchmarkTools.Trial:
  memory estimate:  64.67 MiB
  allocs estimate:  1393408
  --------------
  minimum time:     118.280 ms (10.98% GC)
  median time:      136.361 ms (18.69% GC)
  mean time:        139.964 ms (20.08% GC)
  maximum time:     252.973 ms (52.04% GC)
  --------------
  samples:          36
  evals/sample:     1

obs: este es un poquito mas lento que los otros porque aqui hago el
overapproximation de Xhatk[bi], mientras que no lo hago en los otros
(que lo hacen en Whatk[bi]). para N grande se deberia ver la diferencia que este
es mejor..

=======================================================
=======================================================
algorithm="no_wrapping_effect"
# the input set is empty, but given explicitly as a single void set
D = ContinuousSystem(A, X0, VoidSet(6))

julia> @time compute(400)
System Construction...
elapsed time: 1.7612e-5 seconds
Time Discretization...
elapsed time: 0.000404754 seconds
Reachable LazySets Computation...
no_wrapping_effect
elapsed time: 0.027025543 seconds
Projections...
elapsed time: 0.051585912 seconds
  0.079559 seconds (572.03 k allocations: 28.164 MiB, 13.79% gc time)
0.051585912

BenchmarkTools.Trial:
  memory estimate:  28.16 MiB
  allocs estimate:  572016
  --------------
  minimum time:     61.910 ms (9.74% GC)
  median time:      70.571 ms (10.41% GC)
  mean time:        71.107 ms (12.26% GC)
  maximum time:     85.053 ms (15.63% GC)
  --------------
  samples:          71
  evals/sample:     1
=======================================================
=======================================================
algorithm="no_wrapping_effect"
# the input set is empty, but given explicitly as a time-varying input
U = [VoidSet(6) for i in 1:N]
D = ContinuousSystem(A, X0, U)

julia> @time compute(400)
System Construction...
elapsed time: 2.0093e-5 seconds
Time Discretization...
elapsed time: 0.108007939 seconds
Reachable LazySets Computation...
no_wrapping_effect
elapsed time: 0.026359262 seconds
Projections...
elapsed time: 0.039701652 seconds
  0.174458 seconds (1.21 M allocations: 54.035 MiB, 10.19% gc time)
0.039701652

BenchmarkTools.Trial:
  memory estimate:  54.03 MiB
  allocs estimate:  1208956
  --------------
  minimum time:     160.396 ms (7.19% GC)
  median time:      169.430 ms (9.08% GC)
  mean time:        168.627 ms (9.18% GC)
  maximum time:     181.301 ms (11.25% GC)
  --------------
  samples:          30
  evals/sample:     1
  
=======================================================
=======================================================
algorithm="no_wrapping_effect"
# input set is non-empty and time-varying
U = [BallInf(ones(6), 0.05 * i) for i in 1:N]
D = ContinuousSystem(A, X0, U)

julia> @time compute(400)
System Construction...
elapsed time: 4.7554e-5 seconds
Time Discretization...
elapsed time: 0.11787731 seconds
Reachable LazySets Computation...
no_wrapping_effect
elapsed time: 2.922246121 seconds
Projections...
elapsed time: 0.05933473 seconds
  3.099869 seconds (13.60 M allocations: 715.949 MiB, 10.03% gc time)
0.05933473

BenchmarkTools.Trial:
  memory estimate:  715.95 MiB
  allocs estimate:  13600762
  --------------
  minimum time:     2.968 s (8.40% GC)
  median time:      3.057 s (10.19% GC)
  mean time:        3.057 s (10.19% GC)
  maximum time:     3.146 s (11.89% GC)
  --------------
  samples:          2
  evals/sample:     1
  
=======================================================
========================================================

reach_explicit_block:

julia> @time compute(1500)

System Construction...
elapsed time: 1.5789e-5 seconds
Time Discretization...
elapsed time: 0.000330447 seconds
Reachable States Computation...
elapsed time: 0.058976726 seconds
Projection...
elapsed time: 0.032561023 seconds
Plotting...
elapsed time: 3.918620477 seconds
  4.540615 seconds (1.29 M allocations: 58.947 MiB, 0.48% gc time)
  
  =======================================================
  ========================================================
  
  
  