# Axis log scales, titles, and number/label display (AxisStyle, CanvasStyle,
# LabelChange). No-ROOT tests check validation and the request protocol; native
# tests read the ROOT state recorded by the ROOTPLOTS_AXIS_DUMP hook just before
# each image is saved.

"""Parse a ROOTPLOTS_AXIS_DUMP file: one entry per canvas, with its pads in drawing order."""
function read_axis_dump(path)
    canvases = []
    stack = Any[]
    for line in readlines(path)
        f = split(line, '\t'; keepempty=true)
        if f[1] == "CANVAS"
            push!(canvases, (name=String(f[2]), opttitle=parse(Int, f[3]), pads=Any[], style=Ref{Any}(nothing)))
        elseif f[1] == "GSTYLE"
            v = parse.(Float64, f[3:end])
            last(canvases).style[] = (name=String(f[2]), margins=Tuple(v[1:4]), label_font=Int(v[5]),
                                      label_size=v[6], title_font=Int(v[7]), title_size=v[8],
                                      tickx=Int(v[9]), ticky=Int(v[10]), optstat=Int(v[11]))
        elseif f[1] == "PAD"
            v = parse.(Float64, f[3:end])
            pad = (name=String(f[2]), logx=v[1] == 1, logy=v[2] == 1, logz=v[3] == 1, gridx=v[4] == 1,
                   gridy=v[5] == 1, tickx=Int(v[6]), ticky=Int(v[7]), margins=Tuple(v[8:11]),
                   place=Tuple(v[12:15]),  # canvas NDC: xlow, ylow, width, height
                   title=Ref{Any}(nothing), axes=Dict{String,Any}(), latex=Any[], palette=Ref{Any}(nothing))
            push!(last(canvases).pads, pad); push!(stack, pad)
        elseif f[1] == "ENDPAD"
            pop!(stack)
        elseif f[1] == "LATEX"
            push!(last(stack).latex, (text=String(f[2]), font=parse(Int, f[3]), size=parse(Float64, f[4]),
                                      x=parse(Float64, f[5]), y=parse(Float64, f[6]), align=parse(Int, f[7]),
                                      ndc=f[8] == "1",
                                      # TLatex box (ascent + descent) in pad NDC, measured by ROOT
                                      width=parse(Float64, f[9]), height=parse(Float64, f[10])))
        elseif f[1] == "TITLE"
            last(stack).title[] = (text=String(f[2]), box=parse.(Float64, f[3:6]), fill=String(f[7]),
                                   fill_style=parse(Int, f[8]), border=parse(Int, f[9]),
                                   text_color=String(f[10]), font=parse(Int, f[11]))
        elseif f[1] == "AXIS"
            n = parse.(Int, f[25:28])
            last(stack).axes[String(f[2])] = (title=String(f[3]), ndiv=parse(Int, f[4]), maxdigits=parse(Int, f[5]),
                no_exponent=f[6] == "1", more_log_labels=f[7] == "1", decimals=f[8] == "1",
                center_title=f[9] == "1", center_labels=f[10] == "1", ticks=String(f[11]),
                time_display=f[12] == "1", time_format=String(f[13]), label_color=String(f[14]),
                title_color=String(f[15]), axis_color=String(f[16]), label_size=parse(Float64, f[17]),
                title_size=parse(Float64, f[18]), title_offset=parse(Float64, f[19]),
                label_offset=parse(Float64, f[20]), tick_length=parse(Float64, f[21]),
                label_font=parse(Int, f[22]), title_font=parse(Int, f[23]), labels_option=String(f[24]),
                nbinlabels=n[1], nbins=n[2], first=n[3], last=n[4],
                bin_labels=Dict{Int,String}(), changes=Any[])
        elseif f[1] == "PALETTE"
            v = parse.(Float64, f[2:7])
            last(stack).palette[] = (x1=v[1], y1=v[2], x2=v[3], y2=v[4], title_left=v[5], title_right=v[6])
        elseif f[1] == "BINLABEL"
            last(stack).axes[String(f[2])].bin_labels[parse(Int, f[3])] = String(f[4])
        elseif f[1] == "MODLABEL"
            push!(last(stack).axes[String(f[2])].changes,
                  (index=parse(Int, f[3]), angle=parse(Float64, f[4]), size=parse(Float64, f[5]),
                   align=parse(Int, f[6]), color=String(f[7]), font=parse(Int, f[8]), text=String(f[9])))
        end
    end
    return canvases
end

"""Run f(dir) with ROOTPLOTS_AXIS_DUMP set and return the parsed canvases."""
function with_axis_dump(f)
    mktempdir() do dir
        path = joinpath(dir, "axis.txt")
        withenv("ROOTPLOTS_AXIS_DUMP" => path) do
            f(dir)
        end
        return isfile(path) ? read_axis_dump(path) : []
    end
end
"""Pads of a dumped canvas that own a frame (have axes)."""
framed(canvas) = filter(p -> !isempty(p.axes), canvas.pads)

@testset "AxisStyle/CanvasStyle validation and protocol (no ROOT)" begin
    @test AxisStyle() == AxisStyle() && hash(AxisStyle()) == hash(AxisStyle())
    @test AxisStyle(bin_labels=["a", "b"]).bin_labels == [1 => "a", 2 => "b"]
    @test AxisStyle(bin_labels=[3 => "c", 1 => "a"]).bin_labels == [1 => "a", 3 => "c"]
    @test AxisStyle(bin_labels=Dict(2 => "b")).bin_labels == [2 => "b"]
    @test AxisStyle(bin_labels=["a"]) == AxisStyle(bin_labels=["a"])
    @test AxisStyle(change_labels=LabelChange(1; size=0)).change_labels == [LabelChange(1; size=0)]
    @test AxisStyle(exponent_offset=(0, 0.1)).exponent_offset === (0.0, 0.1)
    @test LabelChange(-1; angle=30, text="last").angle === 30.0
    for bad in (() -> AxisStyle(bin_labels="abc"), () -> AxisStyle(bin_labels=[1.5 => "x"]),
                () -> AxisStyle(change_labels=[1]), () -> AxisStyle(bin_labels=[:a]))
        @test_throws ArgumentError bad()
    end
    v(a; name=:x) = ROOTPlots.validate_axis(a; name)
    @test v(AxisStyle(no_exponent=true, more_log_labels=true, decimals=true, optimize_divisions=false,
                      divisions=505, center_labels=true, label_color="#112233", title_color="#445566",
                      axis_color="#778899", ticks="+-", time_display=true, time_format="%H:%M",
                      time_offset=0.0, exponent_offset=(-0.05, 0.02), bin_labels=["a"], labels_option="v",
                      change_labels=[LabelChange(2; angle=45, size=0.03, align=12, color="#ff0000",
                                                 font=62, text="two")])) isa AxisStyle
    for bad in (AxisStyle(ticks="x"), AxisStyle(ticks="++"), AxisStyle(labels_option="a"),
                AxisStyle(labels_option="v"), AxisStyle(time_format=""), AxisStyle(time_format="a\tb"),
                AxisStyle(time_format="é"), AxisStyle(time_format=repeat("x", 257)),
                AxisStyle(time_offset=Inf), AxisStyle(exponent_offset=(2.0, 0.0)),
                AxisStyle(exponent_offset=(NaN, 0.0)), AxisStyle(bin_labels=String[]),
                AxisStyle(bin_labels=[1 => "a", 1 => "b"]), AxisStyle(bin_labels=[0 => "a"]),
                AxisStyle(bin_labels=[""]), AxisStyle(bin_labels=["a\tb"]), AxisStyle(bin_labels=["a\0"]),
                AxisStyle(label_color="red"), AxisStyle(title_color="#12345"), AxisStyle(axis_color="#gggggg"),
                AxisStyle(change_labels=[LabelChange(0)]), AxisStyle(change_labels=[LabelChange(1001)]),
                AxisStyle(change_labels=[LabelChange(1; angle=Inf)]), AxisStyle(change_labels=[LabelChange(1; size=-1)]),
                AxisStyle(change_labels=[LabelChange(1; align=44)]), AxisStyle(change_labels=[LabelChange(1; font=5)]),
                AxisStyle(change_labels=[LabelChange(1; color="red")]), AxisStyle(change_labels=[LabelChange(1; text="")]),
                AxisStyle(divisions=0), AxisStyle(max_digits=11))
        @test_throws ArgumentError v(bad)
    end
    @test_throws ArgumentError v(AxisStyle(exponent_offset=(0.0, 0.1)); name=:z)

    c(x) = ROOTPlots.validate_canvas(x)
    @test c(CanvasStyle(grid_x=true, grid_y=false, opt_title=true, title_x=0.5, title_y=0.99, title_w=0.4,
                        title_h=0.06, title_align=23, title_font=62, title_font_size=0.05,
                        title_border_size=1, title_fill_color="#eeeeee", title_style=1001,
                        title_text_color="#000088")) isa CanvasStyle
    for bad in (CanvasStyle(title_x=0.5), CanvasStyle(opt_title=false, title_align=23),
                CanvasStyle(opt_title=true, title_x=1.5), CanvasStyle(opt_title=true, title_h=-0.1),
                CanvasStyle(opt_title=true, title_align=44), CanvasStyle(opt_title=true, title_font=5),
                CanvasStyle(opt_title=true, title_font_size=0.0), CanvasStyle(opt_title=true, title_border_size=101),
                CanvasStyle(opt_title=true, title_style=2000), CanvasStyle(opt_title=true, title_fill_color="red"),
                CanvasStyle(opt_title=true, title_text_color="#1234567"))
        @test_throws ArgumentError c(bad)
    end
    merged = ROOTPlots.merge_canvas(CanvasStyle(left_margin=0.14, bottom_margin=0.13),
                                    CanvasStyle(left_margin=0.2, grid_x=true))
    @test merged.left_margin == 0.2 && merged.bottom_margin == 0.13 && merged.grid_x === true

    # Spec types carry the new styles; the original positional forms remain valid.
    @test Graph2DSpec("t", "x", "y", "z", :heatmap, :viridis, 10, false, false, false, (800, 650)).zaxis == AxisStyle()
    @test Graph3DSpec("t", "x", "y", "z", :surface, :viridis, 30., 30., false, false, false, (900, 700)).canvas == CanvasStyle()
    @test PanelSpec("x", "y", "", "", false, false, false, nothing, nothing, nothing).ratio_yaxis == AxisStyle()
    @test PanelSpec("x"; ratio_yaxis=AxisStyle(divisions=503)).ratio_yaxis.divisions == 503
    @test_throws ArgumentError ROOTPlots.validate_nd_spec(Graph2DSpec(zaxis=AxisStyle(ticks="x")))
    @test_throws ArgumentError ROOTPlots.validate_nd_spec(Graph3DSpec(canvas=CanvasStyle(title_x=0.1)))
    @test_throws ArgumentError ROOTPlots.validate_nd_spec(Graph2DSpec(zaxis=AxisStyle(exponent_offset=(0., 0.))))

    # Protocol: one AXIS row (24 scalars) plus its detail rows; default styles add nothing.
    io = IOBuffer()
    ROOTPlots.write_axis(io, "AXIS", "x", AxisStyle(ticks="-", labels_option="d", bin_labels=["μ"],
        time_format="%d/%m", change_labels=[LabelChange(-1; color="#00ff00", text="end")]))
    rows = split(chomp(String(take!(io))), '\n')
    @test length(split(rows[1], '\t')) == 26 && occursin("\tminus\td\t", rows[1])
    @test rows[2] == "AXTIME\tx\t%d/%m" && rows[3] == "AXBIN\tx\t1\t#mu"
    @test rows[4] == "AXCHANGE\tx\t-1\t-1.0\t-1.0\t-1\t#00ff00\t-1\t1\tend"
    grid = Grid2D([0., 1., 2.], [0., 1., 2.], [1. 2.; 3. 4.])
    mktemp() do path, stream
        close(stream)
        ROOTPlots.write_grid_spec(path, grid, Graph2DSpec())
        @test !occursin("SCENE", read(path, String))
        ROOTPlots.write_grid_spec(path, grid, Graph2DSpec(zaxis=AxisStyle(no_exponent=true)))
        text = read(path, String)
        @test occursin("SCENE\t1\nCANVAS\t", text) && occursin("\nAXIS\tz\t", text)
    end

    # Bin labels need a histogram bin axis with enough bins; checked before ROOT runs.
    mktempdir() do dir
        bad = joinpath(dir, "absent", "plot")
        h = Hist1D([0., 1., 2., 3.], [1., 2., 3.], [1., 2., 3.])
        g = GraphSeries([1., 2.], [1., 2.])
        labels3 = AxisStyle(bin_labels=["a", "b", "c"])
        nr = "/not-root"
        for request in (
            () -> render_root([g], GraphSpec(xaxis=labels3), bad; root_config=nr),
            () -> render_root(h, GraphSpec(xaxis=AxisStyle(bin_labels=[4 => "d"])), bad; draw=:e1, root_config=nr),
            () -> render_root(h, GraphSpec(yaxis=labels3), bad; draw="E1", root_config=nr),
            () -> render_root(ModelCurve(:pol1, [1., 1.]; range=(0., 1.)), GraphSpec(xaxis=labels3), bad; root_config=nr),
            () -> render_root(QQGraph([1., 2., 3.], [2., 3., 4.]), GraphSpec(xaxis=labels3), bad; root_config=nr),
            () -> render_root(PolarGraph([0., 1.], [1., 2.]), GraphSpec(xaxis=AxisStyle(label_size=0.05)), bad; root_config=nr),
            () -> render_root(PieChart([1., 2.]), GraphSpec(canvas=CanvasStyle(opt_title=true)), bad; root_config=nr),
            () -> render_root(grid, Graph2DSpec(xaxis=AxisStyle(bin_labels=["a", "b", "c"])), bad; root_config=nr),
            () -> render_root(grid, Graph2DSpec(zaxis=AxisStyle(bin_labels=["a"])), bad; root_config=nr),
            () -> render_root(Graph3D([0., 1., 0.5], [0., 0.5, 1.], [1., 2., 3.]), Graph3DSpec(draw=:points, xaxis=labels3), bad; root_config=nr),
            () -> render_root(ModelSurface(:pol1_2d, [1., 1., 1.]; xrange=(0., 1.), yrange=(0., 1.)), Graph2DSpec(xaxis=labels3), bad; root_config=nr),
            () -> render_root(PolygonHist2D([([0., 1., .5], [0., 0., 1.])], [1.]), Graph2DSpec(xaxis=labels3), bad; root_config=nr),
            () -> render_root(PlotPanel([g], GraphSpec(xaxis=labels3)), bad; root_config=nr),
            () -> render_root(PlotPanel(h, GraphSpec(xaxis=AxisStyle(bin_labels=[9 => "x"]))), bad; root_config=nr),
            () -> render_root(HistStack([h]), GraphSpec(xaxis=AxisStyle(bin_labels=[5 => "x"])), bad; root_config=nr),
            () -> render_root(PlotFigure([PlotPanel(h, GraphSpec(xaxis=AxisStyle(exponent_offset=(0., 0.1)))),
                                          PlotPanel(h, GraphSpec(xaxis=AxisStyle(exponent_offset=(0., 0.2))))]), bad; root_config=nr),
        )
            @test_throws ArgumentError request()
        end
        fit = ROOTPlots.FitResult(:pol0, (0., 3.), [2.], [0.1], fill(0.01, 1, 1), falses(1), 0, true, 3,
                                  1.0, :chi2, 2, 3, 0.0, "Minuit2", false)
        for spec in (GraphSpec(canvas=CanvasStyle(left_margin=0.2)), GraphSpec(xaxis=labels3),
                     GraphSpec(texts=[PlotText(0.5, 0.5, "x")]))
            @test_throws ArgumentError render_root(h, fit, spec, bad; root_config=nr)
        end
        data = Hist1D([0., 1., 2.], [1., 2.], [1., 2.])
        cmp = comparison(data, [Component("MC", data, "#123456")])
        for spec in (PanelSpec("x"; canvas=CanvasStyle(opt_title=true)), PanelSpec("x"; canvas=CanvasStyle(left_margin=0.2)),
                     PanelSpec("x"; xaxis=AxisStyle(bin_labels=["a", "b", "c"])), PanelSpec("x"; ratio_yaxis=AxisStyle(bin_labels=["a"])),
                     PanelSpec("x"; yaxis=AxisStyle(exponent_offset=(0., 0.1)), ratio_yaxis=AxisStyle(exponent_offset=(0., 0.2))))
            @test_throws ArgumentError render_root([cmp], [spec], bad; root_config=nr)
        end
        @test !ispath(dirname(bad))
    end
end

@testset "ROOT axis/pad state reaches ROOT" begin
    if TEST_NATIVE_ROOT
        h = Hist1D(collect(0.:1.:5.), [3., 7., 9., 6., 2.], [3., 7., 9., 6., 2.])
        positive = GraphSeries([1., 10., 100., 1000.], [2., 20., 2e4, 3e5]; label="data")
        labels5 = ["a", "b", "c", "d", "e"]
        number_axis = AxisStyle(no_exponent=true, more_log_labels=true, decimals=true, divisions=505,
                                optimize_divisions=false, center_labels=true, ticks="+-",
                                label_color="#cc3311", title_color="#0044aa", axis_color="#008800",
                                change_labels=[LabelChange(1; angle=30, size=0.05, align=12, color="#ff00ff", font=62, text="first"),
                                               LabelChange(-1; size=0)])
        canvases = with_axis_dump() do dir
            # 1: graph, log x/y, number formats, colours, ticks, grids
            render_root([positive], GraphSpec(logx=true, logy=true, xaxis=number_axis,
                yaxis=AxisStyle(max_digits=3, divisions=510), canvas=CanvasStyle(grid_x=true, grid_y=true)),
                joinpath(dir, "graph"))
            # 2: Hist1D preset (owned DrawFrame frame takes the histogram binning) + vertical bin labels
            render_root(h, GraphSpec(xaxis=AxisStyle(bin_labels=labels5, labels_option="v")), joinpath(dir, "labels"); draw=:e1)
            # 3: default :steps switches to a native TH1D for bin labels; ROOT title with axis titles
            render_root(h, GraphSpec(title="Counts;channel;entries", xaxis=AxisStyle(bin_labels=labels5),
                canvas=CanvasStyle(opt_title=true, title_x=0.3, title_align=13, title_border_size=1,
                                   title_fill_color="#eeeeee", title_style=1001, title_text_color="#000088", title_font=62)),
                joinpath(dir, "steps"))
            # 4: time axis on a raw TH1 draw
            render_root(h, GraphSpec(xaxis=AxisStyle(time_display=true, time_format="%H:%M", time_offset=0)),
                joinpath(dir, "time"); draw="HIST")
            # 5: default plot: no ROOT title, no grids
            render_root(h, GraphSpec(title="Plain"), joinpath(dir, "plain"); draw=:e1)
        end
        @test length(canvases) == 5 && all(c -> c.opttitle == 0, canvases) # gStyle OptTitle is restored
        pad = only(framed(canvases[1])); x = pad.axes["x"]; y = pad.axes["y"]
        @test pad.logx && pad.logy && pad.gridx && pad.gridy
        @test x.no_exponent && x.more_log_labels && x.decimals && x.center_labels && x.ticks == "+-"
        @test x.ndiv == -505 # TAttAxis::SetNdivisions(505, kFALSE)
        @test x.label_color == "#cc3311" && x.title_color == "#0044aa" && x.axis_color == "#008800"
        m = x.changes[1]
        @test (m.index, m.angle, m.align, m.color, m.font, m.text) == (1, 30.0, 12, "#ff00ff", 62, "first")
        @test isapprox(m.size, 0.05; atol=1e-6) # ROOT stores a Float_t
        @test x.changes[2].index == -1 && x.changes[2].size == 0
        @test y.maxdigits == 3 && y.ndiv == 3000510 && !y.no_exponent
        pad = only(framed(canvases[2])); x = pad.axes["x"]
        @test x.nbins == 5 && x.bin_labels == Dict(i => labels5[i] for i in 1:5) && x.labels_option == "v"
        pad = only(framed(canvases[3]))
        @test pad.axes["x"].nbinlabels == 5
        t = pad.title[]
        @test t.text == "Counts" && pad.axes["x"].title == "channel" && pad.axes["y"].title == "entries"
        @test isapprox(t.box[1], 0.3; atol=1e-6) && t.border == 1 && t.fill == "#eeeeee" && t.fill_style == 1001
        @test t.text_color == "#000088" && t.font == 62
        x = only(framed(canvases[4])).axes["x"]
        @test x.time_display && startswith(x.time_format, "%H:%M%F1970-01-01 00:00:00")
        pad = only(framed(canvases[5]))
        @test pad.title[] === nothing && !pad.gridx && !pad.logx && pad.axes["x"].ndiv == 510

        grid = Grid2D([1., 2., 3., 4.], [1., 2., 3., 4.], [1. 20. 300.; 2. 40. 600.; 3. 60. 900.])
        volume = Grid3D([0., 1., 2.], [0., 1., 2.], [1., 2., 3.], reshape(collect(1.:8.), 2, 2, 2))
        h2 = Hist2D([0., 1., 2., 3.], [0., 1., 2.], [1. 2.; 3. 4.; 5. 6.], ones(3, 2))
        zstyle = AxisStyle(no_exponent=true, label_color="#0000ff", title_offset=1.6, more_log_labels=true)
        canvases = with_axis_dump() do dir
            render_root(grid, Graph2DSpec(title="TH2", logz=true, logx=true, zaxis=zstyle,
                canvas=CanvasStyle(opt_title=true, grid_y=true)), joinpath(dir, "colz"))
            render_root(grid, Graph2DSpec(), joinpath(dir, "default2d"))
            render_root(h2, Graph2DSpec(xaxis=AxisStyle(bin_labels=["u", "v", "w"]), zaxis=AxisStyle(decimals=true),
                canvas=CanvasStyle(right_margin=0.25)), joinpath(dir, "h2"); draw="COLZ")
            render_root(volume, Graph3DSpec(draw="BOX2Z", zaxis=AxisStyle(divisions=503), logz=true,
                xaxis=AxisStyle(bin_labels=["p", "q"])), joinpath(dir, "volume"))
            render_root(ModelSurface(:gaus2d, [5., 0., 1., 0., 1.]; xrange=(-2., 2.), yrange=(-2., 2.)),
                Graph2DSpec(logz=true, zaxis=AxisStyle(more_log_labels=true), xaxis=AxisStyle(ticks="-")), joinpath(dir, "tf2"))
            render_root(ModelVolume(:gaus3d, [1., 0., 1., 0., 1., 0., 1.]; xrange=(-2., 2.), yrange=(-2., 2.), zrange=(-2., 2.), level=0.3),
                Graph3DSpec(draw=:iso, xaxis=AxisStyle(divisions=504), canvas=CanvasStyle(opt_title=true), title="TF3"), joinpath(dir, "tf3"))
            render_root(Graph3D([0., 1., 0., 1.], [0., 0., 1., 1.], [1., 2., 3., 4.]),
                Graph3DSpec(draw=:points, yaxis=AxisStyle(label_color="#aa0000"), title="XYZ", canvas=CanvasStyle(opt_title=true)),
                joinpath(dir, "xyz"))
            render_root([grid, grid], Graph2DSpec(zaxis=AxisStyle(no_exponent=true)), joinpath(dir, "overlay");
                draw=["COLZ", "CONT3 SAME"])
        end
        @test length(canvases) == 8
        pad = only(framed(canvases[1]))
        @test pad.logz && pad.logx && pad.gridy && pad.title[].text == "TH2"
        z = pad.axes["z"]
        @test z.no_exponent && z.more_log_labels && z.label_color == "#0000ff" && isapprox(z.title_offset, 1.6; atol=1e-6)
        pad = only(framed(canvases[2]))
        @test pad.title[] === nothing && !pad.logz && isapprox(pad.margins[2], 0.2; atol=1e-6) && isapprox(pad.axes["z"].title_offset, 1.3; atol=1e-6)
        pad = only(framed(canvases[3]))
        @test isapprox(pad.margins[2], 0.25; atol=1e-6) && pad.axes["x"].bin_labels == Dict(1 => "u", 2 => "v", 3 => "w")
        @test pad.axes["z"].decimals
        pad = only(framed(canvases[4]))
        @test pad.logz && pad.axes["z"].ndiv == 503 && pad.axes["x"].nbinlabels == 2
        pad = only(framed(canvases[5]))
        @test pad.logz && pad.axes["z"].more_log_labels && pad.axes["x"].ticks == "-"
        pad = only(framed(canvases[6]))
        @test pad.axes["x"].ndiv == 504 && pad.title[] !== nothing && pad.title[].text == "TF3"
        pad = only(framed(canvases[7]))
        @test pad.axes["y"].label_color == "#aa0000" && pad.title[] !== nothing && pad.title[].text == "XYZ"
        @test only(framed(canvases[8])).axes["z"].no_exponent

        a = Hist1D([1., 2., 3., 4.], [40., 600., 3000.], [40., 600., 3000.])
        stack = HistStack([a, a]; labels=["A", "B"], fill_colors=["#cc6677", "#4477aa"])
        profile = Profile1D([0., 1., 2., 3.], [.2, .4, 1.2, 1.4, 2.5], [1., 3., 2., 5., 4.])
        profile2 = Profile2D([0., 1., 2.], [0., 1., 2.], [.2, .7, 1.4, 1.6], [.2, .5, .3, 1.6], [1., 2., 3., 4.])
        poly = PolygonHist2D([([1., 2., 1.5], [1., 1., 2.])], [5.])
        canvases = with_axis_dump() do dir
            render_root(stack, GraphSpec(logy=true, yaxis=AxisStyle(more_log_labels=true), title="Stack",
                canvas=CanvasStyle(opt_title=true)), joinpath(dir, "stack"))
            render_root(profile, GraphSpec(xaxis=AxisStyle(bin_labels=["lo", "mid", "hi"])), joinpath(dir, "profile"))
            render_root(profile2, Graph2DSpec(zaxis=AxisStyle(label_color="#00aa00")), joinpath(dir, "profile2"))
            render_root(poly, Graph2DSpec(logx=true, logy=true, xaxis=AxisStyle(axis_color="#ff0000")), joinpath(dir, "poly"))
            render_root(QQGraph([1., 2., 3., 4.], [2., 3., 5., 9.]), GraphSpec(canvas=CanvasStyle(grid_x=true),
                yaxis=AxisStyle(divisions=504)), joinpath(dir, "qq"))
            render_root(ModelCurve(:exp, [1., 0.5]; range=(1., 10.)), GraphSpec(logy=true, yaxis=AxisStyle(more_log_labels=true),
                title="TF1", canvas=CanvasStyle(opt_title=true)), joinpath(dir, "tf1"))
            render_root([a, a], GraphSpec(logy=true, yaxis=AxisStyle(no_exponent=true)), joinpath(dir, "overlay");
                draw=["HIST", "E1 SAME"])
            render_root(a, [GraphSeries([1.5, 2.5], [50., 700.])], GraphSpec(xaxis=AxisStyle(bin_labels=["x", "y", "z"]),
                canvas=CanvasStyle(opt_title=true), title="Mixed"), joinpath(dir, "mixed"); draw="HIST")
        end
        @test length(canvases) == 8
        pad = only(framed(canvases[1]))
        @test pad.logy && pad.axes["y"].more_log_labels && pad.title[] !== nothing && pad.title[].text == "Stack"
        @test only(framed(canvases[2])).axes["x"].bin_labels == Dict(1 => "lo", 2 => "mid", 3 => "hi")
        @test only(framed(canvases[3])).axes["z"].label_color == "#00aa00"
        pad = only(framed(canvases[4]))
        @test pad.logx && pad.logy && pad.axes["x"].axis_color == "#ff0000"
        pad = only(framed(canvases[5]))
        @test pad.gridx && pad.axes["y"].ndiv == 504
        pad = only(framed(canvases[6]))
        @test pad.logy && pad.axes["y"].more_log_labels && pad.title[] !== nothing && pad.title[].text == "TF1"
        pad = only(framed(canvases[7]))
        @test pad.logy && pad.axes["y"].no_exponent
        pad = only(framed(canvases[8]))
        @test pad.axes["x"].nbinlabels == 3 && pad.title[] !== nothing && pad.title[].text == "Mixed"

        # Figures: per-panel settings stay per pad (title attributes are restored for the next pad).
        g = GraphSeries([1., 2., 3.], [10., 200., 3000.]; label="g")
        panels = [PlotPanel(g, GraphSpec(logy=true, title="A", yaxis=AxisStyle(more_log_labels=true),
                                         canvas=CanvasStyle(opt_title=true, title_fill_color="#ffeeee", title_style=1001))),
                  PlotPanel(h, GraphSpec(title="B", xaxis=AxisStyle(bin_labels=labels5, labels_option="u"),
                                         canvas=CanvasStyle(grid_y=true))),
                  PlotPanel(h2, Graph2DSpec(title="C", logz=true, zaxis=AxisStyle(no_exponent=true), canvas=CanvasStyle(opt_title=true))),
                  PlotFigure([PlotPanel(h, GraphSpec(xaxis=AxisStyle(divisions=503, optimize_divisions=false))),
                              PlotPanel(grid, Graph2DSpec(logx=true, canvas=CanvasStyle(grid_x=true)))]; layout=(1, 2))]
        canvases = with_axis_dump() do dir
            render_root(PlotFigure(panels; layout=(2, 2), size=(900, 700)), joinpath(dir, "figure"))
        end
        pads = framed(only(canvases))
        @test length(pads) == 5
        @test pads[1].logy && pads[1].axes["y"].more_log_labels && pads[1].title[].text == "A" && pads[1].title[].fill == "#ffeeee"
        @test pads[2].gridy && pads[2].title[] === nothing && pads[2].axes["x"].labels_option == "u" && pads[2].axes["x"].nbins == 5
        @test pads[3].logz && pads[3].axes["z"].no_exponent && pads[3].title[].text == "C" && pads[3].title[].fill == "#ffffff"
        @test pads[4].axes["x"].ndiv == -503 && !pads[4].gridx
        @test pads[5].logx && pads[5].gridx && !pads[5].logz

        # Two-pad views: fit/pull and comparison.
        fh = Hist1D([1., 2., 3., 4., 5.], [4., 6., 5., 7.], [4., 6., 5., 7.])
        fit = fit_root(fh, :pol0; initial=[5.])
        data = Hist1D([1., 2., 3., 4.], [70., 800., 9000.], [70., 800., 9000.])
        mc = Hist1D([1., 2., 3., 4.], [60., 700., 9500.], [60., 700., 9500.])
        cmp = comparison(data, [Component("MC", mc, "#4477aa")])
        canvases = with_axis_dump() do dir
            render_root(fh, fit, GraphSpec(title="Fit", logx=true, xaxis=AxisStyle(divisions=504, label_color="#aa0000"),
                yaxis=AxisStyle(no_exponent=true), canvas=CanvasStyle(grid_x=true, opt_title=true, ticks_y=false)),
                joinpath(dir, "fitpull"))
            render_root([cmp], [PanelSpec("x"; logy=true, xaxis=AxisStyle(bin_labels=["a", "b", "c"]),
                yaxis=AxisStyle(more_log_labels=true), ratio_yaxis=AxisStyle(divisions=505, decimals=true),
                canvas=CanvasStyle(grid_y=true, ticks_x=false))], joinpath(dir, "ratio"))
            render_root([cmp], [PanelSpec("x")], joinpath(dir, "ratio-default"))
        end
        top, bottom = framed(canvases[1])
        @test top.logx && bottom.logx && top.gridx && bottom.gridx && top.ticky == 0
        @test top.title[] !== nothing && top.title[].text == "Fit" && bottom.title[] === nothing
        @test top.axes["x"].ndiv == 504 && top.axes["x"].label_size == 0 && bottom.axes["x"].ndiv == 504
        @test bottom.axes["x"].label_color == "#aa0000" && top.axes["y"].no_exponent
        upper, lower = framed(canvases[2])
        @test upper.logy && upper.gridy && lower.gridy && upper.tickx == 0
        @test upper.axes["y"].more_log_labels && lower.axes["y"].ndiv == 505 && lower.axes["y"].decimals
        @test lower.axes["x"].bin_labels == Dict(1 => "a", 2 => "b", 3 => "c") && upper.axes["x"].label_size == 0
        upper, lower = framed(canvases[3])
        @test !upper.gridy && upper.tickx == 1 && lower.axes["y"].ndiv == 503 && lower.axes["x"].nbinlabels == 0

        # The native parsers reject malformed axis records (exit code 2 = validation error).
        mktempdir() do dir
            exe = ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            input = joinpath(dir, "graph.tsv")
            spec = GraphSpec(xaxis=AxisStyle(bin_labels=["a"], change_labels=[LabelChange(1; text="t")]),
                             canvas=CanvasStyle(grid_x=true))
            ROOTPlots.write_hist_draw_spec(input, h, spec, ROOTPlots.raw_hist_draw_limits(h, spec, GraphStyle(), "H"),
                                           "HIST", "H", GraphStyle())
            valid = read(input, String)
            @test success(pipeline(Cmd([exe, input, joinpath(dir, "ok")]); stdout=devnull, stderr=devnull))
            broken_requests = (replace(valid, "AXBIN\tx\t1\ta" => "AXBIN\tx\t9\ta"),
                               replace(valid, "AXBIN\tx\t1\ta" => "AXBIN\tx\t0\ta"),
                               replace(valid, "AXCHANGE\tx\t1\t" => "AXCHANGE\tx\t0\t"),
                               replace(valid, r"CANVAS(\t[^\t\n]*){9}\t1" => s -> s[1:end-1] * "7"))
            for broken in broken_requests
                @test broken != valid
                write(input, broken)
                result = run(pipeline(ignorestatus(Cmd([exe, input, joinpath(dir, "bad")])); stdout=devnull, stderr=devnull))
                @test result.exitcode == 2
            end
        end
    else
        @test_skip false
    end
end
