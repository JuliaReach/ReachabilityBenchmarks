using ReachabilityBenchmarks, Test
import Aqua

import Pkg
@static if VERSION >= v"1.6"  # TODO make explicit test requirement
    Pkg.add("ExplicitImports")
    import ExplicitImports

    @testset "ExplicitImports tests" begin
        @test isnothing(ExplicitImports.check_all_explicit_imports_are_public(ReachabilityBenchmarks))
        @test isnothing(ExplicitImports.check_all_explicit_imports_via_owners(ReachabilityBenchmarks))
        @test isnothing(ExplicitImports.check_all_qualified_accesses_are_public(ReachabilityBenchmarks))
        @test isnothing(ExplicitImports.check_all_qualified_accesses_via_owners(ReachabilityBenchmarks))
        @test isnothing(ExplicitImports.check_no_implicit_imports(ReachabilityBenchmarks))
        @test isnothing(ExplicitImports.check_no_self_qualified_accesses(ReachabilityBenchmarks))
        @test isnothing(ExplicitImports.check_no_stale_explicit_imports(ReachabilityBenchmarks))
    end
end

@static if VERSION >= v"1.10"
    # JET v0.9.0 (earliest supported version) requires Julia v1.10
    Pkg.add("JET")
    import JET

    @testset "JET tests" begin
        JET.test_package(ReachabilityBenchmarks)
    end
end

@testset "Aqua tests" begin
    Aqua.test_all(ReachabilityBenchmarks)
end
