# Log axes: ROOT clips error bars reaching zero; the fit/pull upper pad is a native TH1D.
using ROOTPlots: graph_limits, fit_pull_limits

"""Primitives (class, name, draw option) of the dumped pads named like pad_name."""
function dumped_primitives(path, pad_name)
    result = Tuple{String,String,String}[]
    current = ""
    for line in readlines(path)
        f = split(line, '\t'; keepempty=true)
        f[1] == "PAD" && (current = String(f[2]))
        f[1] == "DRAWN" && current == pad_name && push!(result, (String(f[2]), String(f[3]), String(f[4])))
    end
    return result
end

@testset "log axes: error bars crossing zero, invalid ranges (no ROOT)" begin
    g = GraphSeries([1.0, 2.0, 3.0], [1.0, 10.0, 100.0]; yerrors=[2.0, 20.0, 5.0])
    (_, (ylo, yhi)) = graph_limits([g], GraphSpec(logy=true))
    @test 0 < ylo < 1.0 && yhi > 105.0     # the range follows the bars' positive parts
    asym = GraphSeries([1.0, 2.0], [5.0, 50.0]; yerrors=(low=[1.0, 60.0], high=[1.0, 1.0]))
    @test graph_limits([asym], GraphSpec(logy=true))[2][1] < 4.0
    gx = GraphSeries([0.5, 2.0], [1.0, 2.0]; xerrors=[1.0, 0.5])
    @test graph_limits([gx], GraphSpec(logx=true))[1][1] > 0
    # Explicit positive limits are kept; nonpositive centres and log ranges stay errors.
    @test graph_limits([g], GraphSpec(logy=true, ylimits=(0.5, 500.0)))[2] == (0.5, 500.0)
    @test_throws ArgumentError graph_limits([GraphSeries([1.0, 2.0], [0.0, 1.0])], GraphSpec(logy=true))
    @test_throws ArgumentError graph_limits([GraphSeries([1.0, 2.0], [-1.0, 1.0])], GraphSpec(logy=true, ylimits=(0.1, 10.0)))
    @test_throws ArgumentError graph_limits([g], GraphSpec(logy=true, ylimits=(0.0, 10.0)))
    @test_throws ArgumentError graph_limits([g], GraphSpec(logy=true, ylimits=(-1.0, 10.0)))
    # Linear axes are unchanged.
    @test graph_limits([g], GraphSpec())[2][1] < -1.0

    # Fit/pull upper pad: empty tails and n=1 bins (error reaching zero) on log y.
    edges = collect(-5.0:1.0:5.0)
    counts = [0.0, 0.0, 1.0, 6.0, 24.0, 38.0, 22.0, 7.0, 1.0, 0.0]
    h = Hist1D(edges, counts, counts)
    data = histogram_points(h)
    curve = GraphSeries(collect(range(-5, 5; length=101)),
                        [38 * exp(-0.5 * x^2) for x in range(-5, 5; length=101)]; draw=:line, errors=:none)
    spec = GraphSpec(logy=true, xlimits=(-5.0, 5.0))
    (_, (lo, hi)) = fit_pull_limits(curve, data, h, spec)
    @test 0 < lo < 1.0                       # n=1 bins (centre 1, bar to 0) are visible
    @test lo > 0.5                           # the Gaussian tail below the data does not set the range
    @test hi > 38 + sqrt(38)
    @test fit_pull_limits(curve, data, h, GraphSpec(logy=true, xlimits=(-5.0, 5.0), ylimits=(0.2, 100.0)))[2] == (0.2, 100.0)
    @test fit_pull_limits(curve, data, h, GraphSpec(xlimits=(-5.0, 5.0)))[2][1] < 0   # linear: unchanged rule
    empty = Hist1D(edges, zeros(10), zeros(10))
    @test_throws ArgumentError fit_pull_limits(curve, histogram_points(empty), empty, spec)
    @test_throws ArgumentError fit_pull_limits(curve, histogram_points(empty), empty,
                                               GraphSpec(logy=true, xlimits=(-5.0, 5.0), ylimits=(0.1, 10.0)))
end

if TEST_NATIVE_ROOT
    @testset "log y fit/pull with empty tails and n=1 bins (native TH1D)" begin
        edges = collect(-5.0:0.5:5.0)
        x = [-2.2; fill(-1.2, 3); fill(-0.7, 9); fill(-0.2, 17); fill(0.3, 16); fill(0.8, 8); fill(1.3, 4); 2.4]
        h = histogram(x, edges).hist
        @test count(==(0), h.sumw) >= 8 && count(==(1), h.sumw) == 2
        fit = fit_root(h, :gaus; initial=[15.0, 0.0, 0.8])
        @test fit.status == 0 && fit.valid
        mktempdir() do dir
            dump = joinpath(dir, "axis.txt")
            withenv("ROOTPLOTS_AXIS_DUMP" => dump) do
                @test isfile(render_root(h, fit, GraphSpec(logy=true, ylabel="Entries"), joinpath(dir, "log")).png)
                @test isfile(render_root(h, fit, GraphSpec(logy=true, ylimits=(0.3, 100.0)), joinpath(dir, "logfixed")).png)
                @test isfile(render_root(h, fit, GraphSpec(), joinpath(dir, "linear")).png)
            end
            upper = dumped_primitives(dump, "fit_top")
            # One TH1D per render (3), drawn by ROOT's histogram painter with E1 X0.
            data = filter(p -> p[2] == "rootplots_fit_data", upper)
            @test length(data) == 3 && all(p -> p[1] == "TH1D" && p[3] == "E1 X0 SAME", data)
            @test !any(p -> p[1] == "TGraphErrors", upper)
            # Invalid log requests still fail in Julia, before ROOT.
            empty = Hist1D(edges, zeros(length(h.sumw)), zeros(length(h.sumw)))
            @test_throws ArgumentError render_root(empty, fit, GraphSpec(logy=true), joinpath(dir, "empty"))
            @test_throws ArgumentError render_root(h, fit, GraphSpec(logy=true, ylimits=(0.0, 10.0)), joinpath(dir, "bad"))
            @test !isfile(joinpath(dir, "empty.png")) && !isfile(joinpath(dir, "bad.png"))
        end
    end

    @testset "log y graphs: error bars crossing zero are drawn and clipped by ROOT" begin
        g = GraphSeries([1.0, 2.0, 3.0], [1.0, 10.0, 100.0]; yerrors=[2.0, 20.0, 5.0], label="data")
        mktempdir() do dir
            @test isfile(render_root([g], GraphSpec(logy=true), joinpath(dir, "graph")).png)
            @test isfile(render_root(PlotPanel([g], GraphSpec(logy=true)), joinpath(dir, "panel")).png)
            @test_throws ArgumentError render_root([GraphSeries([1.0, 2.0], [0.0, 1.0])], GraphSpec(logy=true),
                                                   joinpath(dir, "zero"))
        end
    end
end
