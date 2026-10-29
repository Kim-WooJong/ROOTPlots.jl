# Bent/multi-error graphs, analytic TF2/TF3, non-convex polygons, nested layouts,
# and quiet native ROOT output.

const PNG_MAGIC = UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a]
is_png(path) = isfile(path) && read(path, 8) == PNG_MAGIC
is_pdf(path) = isfile(path) && read(path, 4) == collect(codeunits("%PDF"))
native_fails(exe, path, stem) = !success(pipeline(Cmd([exe, path, stem]); stdout=devnull, stderr=devnull))

@testset "TGraphBentErrors input and serialization" begin
    x = [1., 2., 3.]; y = [2., 3., 2.5]
    bent = BentErrorSeries(x, y; xerrors=[.1, .1, .1], yerrors=(low=[.3, .2, .4], high=[.2, .2, .2]),
                           ybend=(low=[-.1, 0., .1], high=[.1, 0., -.1]), xbend=[.05, .05, .05], label="bent")
    x[1] = 99.
    @test bent.x == [1., 2., 3.] && bent.ybend_low == [-.1, 0., .1] && bent.xbend_high == fill(.05, 3)
    @test bent.eylow == [.3, .2, .4] && bent.exhigh == fill(.1, 3)
    for bad in ((; yerrors=[-.1, .1, .1]), (; ybend=[NaN, 0., 0.]), (; ybend=[0., 0.]),
                (; draw=:fill), (; errors=:unknown), (; xbend=(low=[0., 0., 0.],)), (; label="a\tb"))
        @test_throws ArgumentError BentErrorSeries([1., 2., 3.], [2., 3., 2.5]; bad...)
    end
    @test_throws ArgumentError BentErrorSeries([1., 2.], [1.])
    # Frame ranges include the bent end points, not only the error magnitudes.
    wide = BentErrorSeries([1., 2.], [1., 1.]; yerrors=[.1, .1], ybend=(low=[-.8, 0.], high=[0., .9]))
    xr, _ = ROOTPlots.error_graph_limits(wide, GraphSpec())
    @test xr[1] < .2 && xr[2] > 2.9
    hidden = BentErrorSeries([1., 2.], [1., 1.]; yerrors=[.1, .1], ybend=[5., 5.], errors=:none)
    @test ROOTPlots.error_graph_limits(hidden, GraphSpec())[1][2] < 3
    low = BentErrorSeries([1., 2.], [1., 2.]; xerrors=[.1, .1], xbend=(low=[-1.5, 0.], high=[0., 0.]))
    @test_throws ArgumentError ROOTPlots.error_graph_limits(low, GraphSpec(logy=true))
    io = IOBuffer(); ROOTPlots.write_error_graph(io, bent)
    rows = split(chomp(String(take!(io))), '\n')
    @test length(split(rows[1], '\t')) == 13 && startswith(rows[1], "BENT\t3\tpoints\tbars\tbent")
    @test all(r -> length(split(r, '\t')) == 12 && startswith(r, "BPOINT"), rows[2:end])
end

@testset "TGraphMultiErrors input, sums, and serialization" begin
    x = [1., 2., 3.]; y = [2., 3., 2.5]
    stat = YErrorSet([.3, .4, .5]; label="stat")
    syst = YErrorSet((low=[.4, .3, 0.], high=[.1, .2, .6]); label="syst", display=:boxes, xscale=.5,
                     style=GraphStyle(color="#4477aa", fill_alpha=.4))
    m = MultiErrorSeries(x, y, [stat, syst]; xerrors=[.2, .2, .2], label="data")
    stat.low[1] = 100.
    @test m.sets[1].low == [.3, .4, .5] && length(m.sets) == 2 && m.sets[2].xscale == .5
    @test MultiErrorSeries(x, y, [[.1, .1, .1], (low=[.1, .2, .3], high=[.3, .2, .1])]).sets[2].high == [.3, .2, .1]
    q = total_yerrors(m)
    @test q.low ≈ [.5, .5, .5] && q.high ≈ sqrt.([.3, .4, .5].^2 .+ [.1, .2, .6].^2)
    l = total_yerrors(MultiErrorSeries(x, y, m.sets; sum=:linear))
    @test l.low ≈ [.7, .7, .5] && l.high ≈ [.4, .6, 1.1]
    @test total_yerrors(MultiErrorSeries(x, y, m.sets; sum=:first)).high == [.3, .4, .5]
    @test_throws ArgumentError MultiErrorSeries(x, y, YErrorSet[])
    @test_throws ArgumentError MultiErrorSeries(x, y, fill([.1, .1, .1], 17))
    @test_throws ArgumentError MultiErrorSeries(x, y, [[.1, .1]])
    @test_throws ArgumentError MultiErrorSeries(x, y, [[-.1, .1, .1]])
    @test_throws ArgumentError MultiErrorSeries(x, y, [YErrorSet([.1, .1, .1]; display=:fill)])
    @test_throws ArgumentError MultiErrorSeries(x, y, [YErrorSet([.1, .1, .1]; xscale=0)])
    @test_throws ArgumentError MultiErrorSeries(x, y, [stat]; sum=:mean)
    @test_throws ArgumentError MultiErrorSeries(x, y, [stat]; draw=:bar)
    @test_throws ArgumentError YErrorSet("not errors")
    # Hidden sets do not widen the frame.
    hidden = MultiErrorSeries(x, y, [[.1, .1, .1], YErrorSet(fill(9., 3); display=:none)])
    @test ROOTPlots.error_graph_limits(hidden, GraphSpec())[2][2] < 4
    io = IOBuffer(); ROOTPlots.write_error_graph(io, m)
    rows = split(chomp(String(take!(io))), '\n')
    @test length(split(rows[1], '\t')) == 12 && endswith(rows[1], "\t1")
    @test length(split(rows[2], '\t')) == 11 && occursin("\tboxes\tsyst\t#4477aa", rows[3])
    @test all(r -> length(split(r, '\t')) == 10, rows[4:end])
end

@testset "analytic TF2/TF3 and sampled function grids" begin
    s = ModelSurface(:gaus2d, [2., 0., 1., 1., .5]; xrange=(-2, 2), yrange=(-1, 3))
    @test ROOTPlots.analytic_value(:gaus2d, s.parameters, 0., 1.) == 2.
    @test ROOTPlots.analytic_value(:bigaus, [1., 0., 1., 0., 1., .5], 1., 1.) ≈ exp(-1/1.5)
    @test ROOTPlots.analytic_value(:pol2_3d, collect(1.:10.), 1., 1., 1.) == 55.
    for (model, p) in ((:gaus2d, [1., 0., 0., 0., 1.]), (:bigaus, [1., 0., 1., 0., 1., 1.]),
                       (:pol1_2d, [1., 2.]), (:unknown, [1.]), (:gaus3d, [1., 0., 1., 0., 1., 0., 1.]),
                       (:pol1_2d, [1., NaN, 0.]))
        @test_throws ArgumentError ModelSurface(model, p; xrange=(0, 1), yrange=(0, 1))
    end
    @test_throws ArgumentError ModelSurface(:pol1_2d, [1., 1., 1.]; xrange=(1, 0), yrange=(0, 1))
    @test_throws ArgumentError ModelSurface(:pol1_2d, [1., 1., 1.]; xrange=(0, 1), yrange=(0, 1), npx=2)
    @test_throws ArgumentError ModelVolume(:gaus3d, [1., 0., -1., 0., 1., 0., 1.]; xrange=(0, 1), yrange=(0, 1), zrange=(0, 1))
    @test_throws ArgumentError ModelVolume(:pol1_2d, [1., 1., 1.]; xrange=(0, 1), yrange=(0, 1), zrange=(0, 1))
    @test_throws ArgumentError ModelVolume(:pol1_3d, [1., 1., 1., 1.]; xrange=(0, 1), yrange=(0, 1), zrange=(0, 1), level=Inf)
    @test ROOTPlots.tf2_option(Graph2DSpec()) == "COLZ"
    @test ROOTPlots.tf2_option(Graph2DSpec(draw=:contour)) == "CONT3"
    @test ROOTPlots.tf2_option(Graph3DSpec()) == "SURF2Z"
    @test ROOTPlots.tf2_option(Graph2DSpec(draw="CONT1Z")) == "CONT1Z"
    @test ROOTPlots.tf2_option(Graph3DSpec(draw="LEGO2")) == "LEGO2"
    for bad in (Graph2DSpec(draw="SURF1"), Graph3DSpec(draw="COLZ"), Graph3DSpec(draw="GLSURF"),
                Graph2DSpec(draw="CONT1Z; x"), Graph2DSpec(draw="BOX"), Graph3DSpec(draw=:points))
        @test_throws ArgumentError ROOTPlots.tf2_option(bad)
    end
    v = ModelVolume(:pol2_3d, [-1., 0, 0, 0, 1, 1, 1, 0, 0, 0]; xrange=(-2, 2), yrange=(-2, 2), zrange=(-2, 2))
    @test ROOTPlots.tf3_mode(Graph3DSpec(draw=:iso)) == :iso
    for bad in (Graph3DSpec(), Graph3DSpec(draw="BOX2Z"), Graph3DSpec(draw=:iso, logz=true))
        @test_throws ArgumentError ROOTPlots.validate_analytic_spec(v, bad)
    end
    @test_throws ArgumentError ROOTPlots.validate_analytic_spec(v, Graph2DSpec())
    @test_throws ArgumentError ROOTPlots.validate_analytic_spec(ModelSurface(:pol1_2d, [0., 1., 0.]; xrange=(-1, 1), yrange=(0, 1)),
                                                                Graph2DSpec(logz=true))
    @test_throws ArgumentError ROOTPlots.validate_analytic_spec(s, Graph2DSpec(logx=true))
    mktempdir() do dir
        rows = readlines(ROOTPlots.write_analytic_spec(joinpath(dir, "tf2.tsv"), s, Graph2DSpec()))
        @test length(split(rows[2], '\t')) == 11 && length(split(rows[4], '\t')) == 12 + 5
        rows = readlines(ROOTPlots.write_analytic_spec(joinpath(dir, "tf3.tsv"), v, Graph3DSpec(draw=:box)))
        @test length(split(rows[2], '\t')) == 7 && length(split(rows[4], '\t')) == 16 + 10
    end
    g = function_grid((x, y) -> x + 10y, (0, 2), (0, 1); n=(2, 4))
    @test g isa Grid2D && g.xedges == [0., 1., 2.] && g.values[2, 1] == 1.5 + 10 * .125
    g3 = function_grid((x, y, z) -> x * y * z, (0, 1), (0, 1), (0, 2); n=(1, 1, 2))
    @test g3 isa Grid3D && g3.values[1, 1, 2] == .5 * .5 * 1.5
    @test_throws ArgumentError function_grid((x, y) -> NaN, (0, 1), (0, 1))
    @test_throws ArgumentError function_grid((x, y) -> "a", (0, 1), (0, 1))
    @test_throws ArgumentError function_grid((x, y) -> x, (0, 1), (0, 1); n=(0, 3))
end

@testset "simple (non-convex) polygons" begin
    L = ([0., 2, 2, 1, 1, 0], [0., 0, 1, 1, 3, 3])
    star = ([cos(pi/2 + i*pi/5) * (isodd(i) ? .4 : 1.) for i in 0:9],
            [sin(pi/2 + i*pi/5) * (isodd(i) ? .4 : 1.) for i in 0:9])
    @test ROOTPlots.simple_polygon(L...) == L
    @test length(ROOTPlots.simple_polygon(star...)[1]) == 10
    @test ROOTPlots.simple_polygon([0., 1, 0, 0], [0., 0, 1, 0]) == ([0., 1, 0], [0., 0, 1])
    bad = [([0., 1, 0, 1], [0., 1, 1, 0]),          # bow tie
           ([0., 1, 2], [0., 0, 0]),                # collinear, zero area
           ([0., 1], [0., 1]),                      # too few vertices
           ([0., 1, NaN], [0., 0, 1]),              # non-finite
           ([0., 2, 2, 1, 1, 0], [0., 0, 2, 0, 3, 3]), # vertex touching a non-adjacent edge
           ([0., 2, 1, 3, 0], [0., 0, 0, 1, 1]),    # edge folding back onto the previous one
           ([0., 1, 1, 0], [0., 0, 0, 1])]          # repeated neighbour
    for (x, y) in bad
        @test_throws ArgumentError ROOTPlots.simple_polygon(x, y)
    end
    @test_throws ArgumentError ROOTPlots.simple_polygon(collect(1.:257.), collect(1.:257.))
    poly = PolygonHist2D([L, star], [3., 5.])
    @test size(poly.polygons[2]) == (10, 2)
    @test_throws ArgumentError PolygonHist2D([([0., 1, 0, 1], [0., 1, 1, 0])], [1.])
    p = PlotPolygon(L...; fill_color="#cc3311", width=2)
    record = ROOTPlots.primitive_record(p, GraphSpec())
    @test startswith(record, "PRIMITIVE\tPOLYGON\t6\t#000000\t2\t#cc3311\t1001\t0.5\t0.0\t0.0\t2.0")
    @test length(split(record, '\t')) == 8 + 12
    @test_throws ArgumentError ROOTPlots.primitive_record(p, GraphSpec(logx=true))
    @test_throws ArgumentError PlotPolygon(L...; fill_color="red")
    @test_throws ArgumentError PlotPolygon(L...; fill_alpha=2)
    @test_throws ArgumentError PlotPolygon([0., 1, 0, 1], [0., 1, 1, 0])
end

@testset "nested figure layouts" begin
    h = Hist1D([0., 1., 2.], [2., 4.], [2., 4.])
    panel = PlotPanel(h, GraphSpec())
    inner = PlotFigure([panel, panel]; layout=(2, 1))
    outer = PlotFigure([panel, inner]; layout=(1, 2), widths=[1, 2])
    @test outer.panels[2] isa PlotFigure && outer.widths == [1., 2.] && outer.heights === nothing
    @test length(ROOTPlots.figure_leaves(outer)) == 3 && ROOTPlots.figure_depth(outer) == 2
    @test PlotFigure([panel inner]).layout == (1, 2)
    @test_throws ArgumentError PlotFigure([panel, 1])
    @test_throws ArgumentError PlotFigure([panel]; widths=[1, 2])
    @test_throws ArgumentError PlotFigure([panel]; heights=[0])
    @test_throws ArgumentError PlotFigure([panel, panel]; widths=[1, 1000])
    deep = PlotFigure([panel])
    for _ in 1:3 deep = PlotFigure([deep]) end
    @test ROOTPlots.figure_depth(deep) == 4
    @test_throws ArgumentError PlotFigure([deep])
    sixteen = PlotFigure(fill(panel, 16); layout=(4, 4))
    @test_throws ArgumentError PlotFigure(fill(sixteen, 5); layout=(1, 5))
    payloads = [["P$i"] for i in 1:3]; decorations = [String[] for _ in 1:3]
    io = IOBuffer()
    ROOTPlots.write_figure_grid(io, outer, payloads, decorations, Ref(0))
    rows = split(chomp(String(take!(io))), '\n')
    @test rows == ["WEIGHTS\t1.0\t2.0\t1.0", "PANEL\t1\t1\t0", "P1", "SUBGRID\t2\t2\t1\t2",
                   "PANEL\t1\t1\t0", "P2", "PANEL\t2\t1\t0", "P3", "ENDGRID\t2"]
    # A regular figure keeps the original flat record stream.
    io = IOBuffer()
    ROOTPlots.write_figure_grid(io, PlotFigure([panel, panel]), payloads[1:2], decorations[1:2], Ref(0))
    @test split(chomp(String(take!(io))), '\n') == ["PANEL\t1\t1\t0", "P1", "PANEL\t2\t1\t0", "P2"]
end

if TEST_NATIVE_ROOT
    @testset "native bent/multi errors, TF2/TF3, polygons, nested pads" begin
        exe = ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
        mktempdir() do dir
            x = collect(1.:5.); y = [2., 3., 2.5, 3.5, 3.]
            spec = GraphSpec(xlabel="x", ylabel="y")
            bent = BentErrorSeries(x, y; xerrors=fill(.2, 5), yerrors=(low=fill(.4, 5), high=fill(.3, 5)),
                                   ybend=(low=fill(-.1, 5), high=fill(.1, 5)), label="bent")
            files = render_root(bent, spec, joinpath(dir, "bent"); format=:both)
            @test is_png(files.png) && is_pdf(files.pdf)
            for (i, mode) in enumerate((:band, :boxes, :caps, :nocaps, :arrows, :none))
                b = BentErrorSeries(x, y; yerrors=fill(.3, 5), ybend=fill(.1, 5), errors=mode, draw=:linepoints)
                @test is_png(render_root(b, spec, joinpath(dir, "bent-$i")).png)
            end
            sets = [YErrorSet(fill(.2, 5); label="stat"),
                    YErrorSet(fill(.5, 5); label="syst", display=:boxes, xscale=.5,
                              style=GraphStyle(color="#4477aa", fill_alpha=.4))]
            multi = MultiErrorSeries(x, y, sets; xerrors=fill(.3, 5), label="data")
            files = render_root(multi, spec, joinpath(dir, "multi"); format=:both)
            @test is_png(files.png) && is_pdf(files.pdf)
            for (i, d) in enumerate(ROOTPlots.ERROR_SET_DISPLAYS)
                m = MultiErrorSeries(x, y, [YErrorSet(fill(.3, 5); display=d), [.1, .1, .1, .1, .1]];
                                     draw=:linepoints, sum=:linear)
                @test is_png(render_root(m, spec, joinpath(dir, "multi-$i")).png)
            end
            # Only boxes: ROOT's painter omits markers, so they are redrawn on top.
            boxes = MultiErrorSeries(x, y, [YErrorSet(fill(.3, 5); display=:boxes)]; xerrors=fill(.2, 5))
            @test is_png(render_root(boxes, spec, joinpath(dir, "multi-boxes")).png)
            mixed = render_root([GraphSeries(x, y .- 1; label="plain"), bent, multi,
                                 ModelCurve(:pol1, [1., .4]; range=(1., 5.), label="TF1")],
                                GraphSpec(legend_spec=LegendSpec(entries=[LegendEntry(:series, 3, "multi")])),
                                joinpath(dir, "mixed"); format=:both)
            @test is_png(mixed.png) && is_pdf(mixed.pdf)
            @test is_png(render_root([bent], GraphSpec(logy=true), joinpath(dir, "bent-log")).png)

            s = ModelSurface(:bigaus, [1., 0., 1., 0., .7, .5]; xrange=(-3, 3), yrange=(-3, 3))
            for (i, sp) in enumerate((Graph2DSpec(), Graph2DSpec(draw=:contour, contours=8),
                                      Graph2DSpec(draw="CONT1Z"), Graph3DSpec(), Graph3DSpec(draw="SURF1"),
                                      Graph3DSpec(draw="LEGO2"), Graph2DSpec(draw="COLZ", logz=false, palette=:bird)))
                files = render_root(s, sp, joinpath(dir, "tf2-$i"); format=:both)
                @test is_png(files.png) && is_pdf(files.pdf)
            end
            positive = ModelSurface(:pol1_2d, [1., 1., 1.]; xrange=(1, 10), yrange=(1, 10))
            @test is_png(render_root(positive, Graph2DSpec(logx=true, logy=true, logz=true), joinpath(dir, "tf2-log")).png)
            v = ModelVolume(:gaus3d, [1., 0., 1., 0., .7, 0., 1.3]; xrange=(-3, 3), yrange=(-3, 3), zrange=(-3, 3), level=.3)
            for mode in (:iso, :box)
                files = render_root(v, Graph3DSpec(draw=mode), joinpath(dir, "tf3-$mode"); format=:both)
                @test is_png(files.png) && is_pdf(files.pdf) && filesize(files.png) > 5000
            end
            grid = function_grid((x, y) -> sin(x) * cos(y), (0, pi), (0, pi); n=(30, 30))
            @test is_png(render_root(grid, Graph2DSpec(), joinpath(dir, "function-grid")).png)

            # A non-convex L must not be filled like its convex hull.
            L = ([0., 2, 2, 1, 1, 0], [0., 0, 1, 1, 3, 3])
            hull = ([0., 2, 2, 1, 0], [0., 0, 1, 3, 3])
            small = Graph2DSpec(size=(400, 400))
            a = render_root(PolygonHist2D([L], [1.]), small, joinpath(dir, "L")).png
            b = render_root(PolygonHist2D([hull], [1.]), small, joinpath(dir, "hull")).png
            @test ROOTPlots.image_difference(a, b).changed_fraction > .01
            primitive = PlotPanel([GraphSeries([0., 3.], [0., 3.])], GraphSpec();
                                  primitives=[PlotPolygon(L...), PlotPolygon(L[1] .+ 1, L[2] .+ 1; fill_style=0)])
            files = render_root(primitive, joinpath(dir, "polygon-primitive"); format=:both)
            @test is_png(files.png) && is_pdf(files.pdf)
            path = joinpath(dir, "poly.tsv")
            ROOTPlots.write_advanced_plot(path, PolygonHist2D([L], [1.]), Graph2DSpec())
            write(path, replace(read(path, String), "VERTEX\t3\t2.0\t1.0" => "VERTEX\t3\t0.5\t-1.0"))
            @test native_fails(exe, path, joinpath(dir, "bad-poly"))

            panel = PlotPanel([bent], GraphSpec(title="bent"))
            inner = PlotFigure([PlotPanel(s, Graph2DSpec(title="TF2")), PlotPanel([multi], GraphSpec(title="multi"))];
                               layout=(2, 1), heights=[1, 2])
            nested = PlotFigure([panel inner; PlotPanel(PolygonHist2D([L], [1.]), Graph2DSpec()) PlotPanel(v, Graph3DSpec(draw=:iso))];
                                size=(1000, 800), widths=[2, 3])
            files = render_root(nested, joinpath(dir, "nested"); format=:both)
            @test is_png(files.png) && is_pdf(files.pdf) && filesize(files.png) > 20000
            flat = render_root(PlotFigure([panel, panel]), joinpath(dir, "flat")).png
            @test is_png(flat)
            malformed = joinpath(dir, "nested.tsv")
            write(malformed, "ROOTPlotsGraph\t1\nDRAWLAYOUT\t800\t600\t1\t1\t1\nSUBGRID\t1\t1\t1\t1\nEND\t1\n")
            @test native_fails(exe, malformed, joinpath(dir, "bad-nested"))
        end
    end

    @testset "quiet ROOT output by default" begin
        mktempdir() do dir
            g = GraphSeries([0., 1.], [1., 2.])
            path = joinpath(dir, "graph.tsv")
            ROOTPlots.write_graph_spec(path, [g], GraphSpec(), ROOTPlots.graph_limits([g], GraphSpec()))
            exe = ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            function captured(cmd)
                err = IOBuffer()
                run(pipeline(cmd; stdout=devnull, stderr=err))
                return String(take!(err))
            end
            quiet = withenv("ROOTPLOTS_ROOT_VERBOSE" => nothing) do
                captured(Cmd([exe, path, joinpath(dir, "quiet")]))
            end
            @test !occursin("TCanvas::Print", quiet) && isfile(joinpath(dir, "quiet.png"))
            verbose = captured(addenv(Cmd([exe, path, joinpath(dir, "verbose")]), "ROOTPLOTS_ROOT_VERBOSE" => "1"))
            @test occursin("TCanvas::Print", verbose)
            renderer = build_root().renderer
            h = Hist1D([0., 1., 2.], [2., 3.], [2., 3.])
            spec = ROOTPlots.write_root_spec(joinpath(dir, "control.tsv"), [comparison(h, [Component("M", h, "#336699")])],
                                             [PanelSpec("x")], (.5, 1.5))
            @test !occursin("TCanvas::Print", withenv(() -> captured(Cmd([renderer, spec, joinpath(dir, "c1")])),
                                                      "ROOTPLOTS_ROOT_VERBOSE" => nothing))
            @test occursin("TCanvas::Print", captured(addenv(Cmd([renderer, spec, joinpath(dir, "c2")]),
                                                             "ROOTPLOTS_ROOT_VERBOSE" => "1")))
            # Through the public API, with the child's stderr redirected to a file.
            log = joinpath(dir, "stderr.txt")
            withenv("ROOTPLOTS_ROOT_VERBOSE" => nothing) do
                open(log, "w") do io
                    redirect_stderr(io) do
                        render_root([g], GraphSpec(), joinpath(dir, "api"); format=:both)
                        render_root(PlotPanel(h, GraphSpec()), joinpath(dir, "api-layout"))
                    end
                end
            end
            @test !occursin("TCanvas::Print", read(log, String))
            @test isfile(joinpath(dir, "api.png")) && isfile(joinpath(dir, "api.pdf"))
        end
    end
end
