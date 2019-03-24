using Dates

# print current time
now()

# Van der Pol benchmark
include("VanDerPol/vanderpol_benchmark.jl")

# Laub-Loomis benchmark
include("LaubLoomis/laubloomis_benchmark.jl")

# Quadrotor benchmark
include("Quadrotor/quadrotor_benchmark.jl")

nothing
