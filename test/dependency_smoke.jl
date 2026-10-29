# Run in a resolved package environment, checking its transitive dependencies.
using Test
using Pkg
using ROOTPlots
using CSV, DataFrames

@testset "standalone package without Julia plotting packages" begin
    names = Set(info.name for info in values(Pkg.dependencies()))
    @test isempty(intersect(names, Set(["CairoMakie", "Makie", "Plots"])))
    @test histogram(CSV.read(IOBuffer("x,w\n0.2,2\n1.2,-3\n"), DataFrame),
                    :x, [0., 1., 2.]; weight=:w).hist.sumw == [2., -3.]
    h = histogram([-.2, .2, .7, 1.2, 2.], [0., 1., 2.]; weights=[2., 1., -2., 3., 4.])
    @test h.hist.sumw == [-1., 3.]
    @test h.hist.sumw2 == [5., 9.]
    @test h.underflow == (2., 4.) && h.overflow == (4., 16.)
    nd = histogram(([.2, 1.2], [.4, 1.4]), ([0., 1., 2.], [0., 1., 2.])).hist
    @test projection(nd, :x).sumw == [1., 1.]
    fixture = FitResult(:pol0, (0., 2.), [1.], [.5], fill(.25, 1, 1), falses(1),
                        0, true, 3, 3.841458820694124, :chi2, 1, 2, 0., "fixture", false)
    @test chi2_pvalue(fixture) ≈ .05 atol=1e-10
end
