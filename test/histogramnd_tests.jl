using Test
using DataFrames
using ROOTPlots
isdefined(ROOTPlots, :Hist2D) || Base.include(ROOTPlots, joinpath(@__DIR__, "..", "src", "histogramnd.jl"))

@testset "Owned multidimensional histograms" begin
    xe, ye, ze = [0., 1, 3], [-2., 0, 4, 9], [1., 2, 5, 7, 10]
    w = reshape(collect(1.:6.), 2, 3)
    v = 2 .* w
    h = ROOTPlots.Hist2D(xe, ye, w, v)
    @test size(h.sumw) == (2, 3)
    @test h.xedges == xe && h.yedges == ye
    xe[1] = -10; ye[1] = -20; w[1] = 99; v[1] = 100
    @test h.xedges[1] == 0 && h.yedges[1] == -2
    @test h.sumw[1] == 1 && h.sumw2[1] == 2
    @test_throws ArgumentError ROOTPlots.Hist2D([0, 0], [0, 1], ones(1,1), ones(1,1))
    @test_throws ArgumentError ROOTPlots.Hist2D([0, 1], [0, Inf], ones(1,1), ones(1,1))
    @test_throws ArgumentError ROOTPlots.Hist2D([0, 1], [0, 1], ones(2,1), ones(2,1))
    @test_throws ArgumentError ROOTPlots.Hist2D([0, 1], [0, 1], fill(NaN,1,1), ones(1,1))
    @test_throws ArgumentError ROOTPlots.Hist2D([0, 1], [0, 1], ones(1,1), -ones(1,1))
    a = reshape(collect(1.:24.), 2, 3, 4)
    h3 = ROOTPlots.Hist3D([0,1,3], [-2,0,4,9], ze, a, 3 .* a)
    @test size(h3.sumw) == (2,3,4)
    a[1] = 999; ze[1] = -999
    @test h3.sumw[1] == 1 && h3.zedges[1] == 1
    @test_throws ArgumentError ROOTPlots.Hist3D([0,1], [0,1], [0,1], ones(1,1,2), ones(1,1,2))
    @test_throws ArgumentError ROOTPlots.Hist3D([0,1], [0,1], [0,1], ones(1,1,1), fill(Inf,1,1,1))
    @test ROOTPlots.validate_stats_hist_nd(h) == (2,3)
    @test ROOTPlots.validate_stats_hist_nd(h3) == (2,3,4)
    h.xedges[2] = h.xedges[1]
    @test_throws ArgumentError ROOTPlots.validate_stats_hist_nd(h)
    @test_throws ArgumentError ROOTPlots.projection(h, :x)
    h3.sumw2[1] = -1
    @test_throws ArgumentError ROOTPlots.projection(h3, (:x,:y))
end

@testset "Vector CPU histogram and ROOT boundaries" begin
    x, w, e = [-1.,0,0.5,1,3,4], [2.,1,-3,4,5,-2], [0.,1,3]
    r = histogram(x, e; weights=w, threads=:off)
    df = histogram(DataFrame(x=x,w=w), :x, e; weight=:w, threads=:off)
    @test r.hist.sumw == df.hist.sumw == [-2.,4]
    @test r.hist.sumw2 == df.hist.sumw2 == [10.,16]
    @test r.underflow == df.underflow == (2.,4.)
    @test r.overflow == df.overflow == (3.,29.)
    @test histogram(Float64[], [0,1]).hist.sumw == [0.]
    @test_throws ArgumentError histogram(x, e; weights=[1.])
    @test_throws ArgumentError histogram(x, e; threads=:bad)
    @test_throws ArgumentError histogram([missing], e)
    @test_throws ArgumentError histogram([NaN], e)
    @test_throws ArgumentError histogram([0.], e; weights=[Inf])
    @test_throws ArgumentError histogram([0.], e; weights=[1e200])
end

@testset "2D/3D flow bins and weighted filling" begin
    xe, ye, ze = [0.,1,3], [-2.,0,4,9], [1.,2,5]
    x = [-1.,0,0.5,1,3,0,0,3]
    y = [-3.,-2,-1,4,0,-3,9,9]
    w = [2.,1,-3,4,5,-2,7,8]
    r = histogram((x,y), (xe,ye); weights=w)
    expected = zeros(4,5); expectedvar = zeros(4,5)
    slots = [(1,1),(2,2),(2,2),(3,4),(4,3),(2,1),(2,5),(4,5)]
    for (slot, weight) in zip(slots,w)
        expected[slot...] += weight
        expectedvar[slot...] += weight^2
    end
    @test size(r.hist.sumw) == (2,3)
    @test r.hist.sumw == expected[2:3,2:4]
    @test r.hist.sumw2 == expectedvar[2:3,2:4]
    @test r.bins_with_flow.sumw == expected
    @test r.bins_with_flow.sumw2 == expectedvar
    r.hist.sumw[1] = 100
    @test r.bins_with_flow.sumw[2,2] == -2
    z = [0.,1,2,5,3,1,2,5]
    r3 = histogram((x,y,z), (xe,ye,ze); weights=w)
    expected3 = zeros(4,5,4); expectedvar3 = zeros(4,5,4)
    zslots = [1,2,3,4,3,2,3,4]
    for (slot, zs, weight) in zip(slots,zslots,w)
        expected3[slot...,zs] += weight
        expectedvar3[slot...,zs] += weight^2
    end
    @test size(r3.hist.sumw) == (2,3,2)
    @test r3.hist.sumw == expected3[2:3,2:4,2:3]
    @test r3.hist.sumw2 == expectedvar3[2:3,2:4,2:3]
    @test r3.bins_with_flow.sumw == expected3
    @test r3.bins_with_flow.sumw2 == expectedvar3
    df = DataFrame(x=x,y=y,z=z,w=w)
    @test histogram(df, (:x,:y), (xe,ye); weight=:w).bins_with_flow == r.bins_with_flow
    @test histogram(df, (:x,:y,:z), (xe,ye,ze); weight=:w).bins_with_flow == r3.bins_with_flow
    @test sum(histogram((x,y), (xe,ye)).bins_with_flow.sumw) == length(x)
    @test size(histogram((Float64[],Float64[]), (xe,ye)).hist.sumw) == (2,3)
    @test_throws ArgumentError histogram((x,y[1:2]), (xe,ye))
    @test_throws ArgumentError histogram((x,y), (xe,ye); weights=[1.])
    @test_throws ArgumentError histogram((x,y), (xe,ye); threads=:bad)
    @test_throws ArgumentError histogram(([0.],[NaN]), (xe,ye))
    @test_throws ArgumentError histogram(([0.],[0.]), (xe,ye); weights=[Inf])
    @test_throws ArgumentError histogram(([0.],[0.]), (xe,ye); weights=[1e200])
    @test_throws ArgumentError histogram((fill(0.,10),fill(0.,10)), (xe,ye); weights=fill(1e154,10))
    @test_throws ArgumentError histogram(df, (:x,:absent), (xe,ye))
    @test_throws ArgumentError histogram(df, (:x,:y), (xe,ye); weight=:absent)
end

@testset "Visible-bin projections preserve requested orientation" begin
    a = reshape(collect(1.:24.), 2,3,4)
    h = ROOTPlots.Hist3D([0,1,3], [-2,0,4,9], [1,2,5,7,10], a, 2 .* a)
    edges = (h.xedges,h.yedges,h.zedges)
    for (i, axis) in enumerate((:x,:y,:z))
        p = ROOTPlots.projection(h, axis)
        dims = Tuple(j for j in 1:3 if j != i)
        @test p.edges == edges[i]
        @test p.sumw == vec(sum(a; dims))
        @test p.sumw2 == vec(sum(2 .* a; dims))
    end
    for (i, j) in ((1,2),(2,1),(1,3),(3,1),(2,3),(3,2))
        axes = ((:x,:y,:z)[i], (:x,:y,:z)[j])
        p = ROOTPlots.projection(h, axes)
        k = only(setdiff(1:3,(i,j)))
        expected = dropdims(sum(a; dims=k); dims=k)
        i > j && (expected = permutedims(expected))
        @test p.xedges == edges[i] && p.yedges == edges[j]
        @test p.sumw == expected
        @test p.sumw2 == 2 .* expected
    end
    h2 = ROOTPlots.Hist2D([0,1,3], [-2,0,4,9], reshape(1.:6.,2,3), reshape(2.:2.:12.,2,3))
    @test ROOTPlots.projection(h2,:x).sumw == [9.,12.]
    @test ROOTPlots.projection(h2,:y).sumw == [3.,7.,11.]
    @test ROOTPlots.projection(h2,:y).sumw2 == [6.,14.,22.]
    p = ROOTPlots.projection(h,:x)
    p.edges[1] = -99; p.sumw[1] = -99
    @test h.xedges[1] == 0 && h.sumw[1] == 1
    @test_throws ArgumentError ROOTPlots.projection(h,:q)
    @test_throws ArgumentError ROOTPlots.projection(h2,:z)
    @test_throws ArgumentError ROOTPlots.projection(h,(:x,:x))
end

@testset "Deterministic chunked filling and bounded private scratch" begin
    n = 3 * ROOTPlots.HIST_CHUNK_SIZE + 17
    x = [mod(i,17)/5 - 1 for i in 1:n]
    y = [mod(i,13)/3 - 1 for i in 1:n]
    z = [mod(i,11)/2 - 1 for i in 1:n]
    weights = [(-1.)^i * (1 + mod(i,7)/10) for i in 1:n]
    edges = ([-1.,0,1,3], [-1.,0,2,4], [-1.,0,2,5])
    for dim in (2,3)
        values, bins = dim == 2 ? (x,y) : (x,y,z), edges[1:dim]
        off = histogram(values,bins; weights, threads=:off)
        on = histogram(values,bins; weights, threads=:on)
        @test isequal(off.hist.sumw,on.hist.sumw)
        @test isequal(off.hist.sumw2,on.hist.sumw2)
        @test isequal(off.bins_with_flow,on.bins_with_flow)
    end
    one_off = histogram(x,edges[1]; weights,threads=:off)
    one_on = histogram(x,edges[1]; weights,threads=:on)
    @test isequal(one_off.hist.sumw,one_on.hist.sumw)
    @test isequal(one_off.hist.sumw2,one_on.hist.sumw2)
    @test isequal(one_off.underflow,one_on.underflow) && isequal(one_off.overflow,one_on.overflow)
    overflowweights = zeros(n)
    overflowweights[1] = overflowweights[end] = 1e154
    @test_throws ArgumentError histogram((fill(0.,n),fill(0.,n)),([0,1],[0,1]); weights=overflowweights)
    # At this shape, four ordinary chunks exceed 64 MiB; choose only three.
    shape = (1100,1100)
    chunks = ROOTPlots.hist_nd_chunk_count(n,shape)
    @test chunks == 3
    @test chunks * 16 * prod(shape) <= ROOTPlots.HIST_SCRATCH_BYTES
    bigedges = collect(0.:1098.)
    smallx, smally = fill(0.5,n), fill(1.5,n)
    off = histogram((smallx,smally),(bigedges,bigedges); threads=:off)
    on = histogram((smallx,smally),(bigedges,bigedges); threads=:on)
    @test isequal(off.bins_with_flow,on.bins_with_flow)
    @test off.hist.sumw[1,2] == n
    @test ROOTPlots.hist_nd_chunk_count(n,(10^9,10^9,10^9)) == 1
end
