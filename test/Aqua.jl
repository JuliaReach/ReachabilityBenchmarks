using ReachabilityBenchmarks, Test
import Aqua

@testset "Aqua tests" begin
    Aqua.test_all(ReachabilityBenchmarks)
end
