using Test, ReachabilityBenchmarks

@testset "Relative path macro" begin
    file = @relpath "my_data.dat"
    @test file == joinpath(@__DIR__, "my_data.dat")
end

include("Aqua.jl")
