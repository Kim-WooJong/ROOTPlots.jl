# ROOT styles (gROOT->SetStyle) and CMS/ATLAS labels (ROOT TLatex). Native
# tests read the ROOT state recorded by the ROOTPLOTS_AXIS_DUMP hook (see
# axis_format_tests.jl): the selected gStyle, pad margins, axis attributes, and
# every TLatex in each pad.

scene_rows(spec) = (io = IOBuffer(); ROOTPlots.write_scene(io, spec); split(String(take!(io)), '\n'; keepempty=false))
label_row(spec) = (rows = filter(r -> startswith(r, "EXPLABEL"), scene_rows(spec)); isempty(rows) ? nothing : split(only(rows), '\t'))
latex(pad, text) = filter(t -> t.text == text, pad.latex)
# ROOT stores sizes and margins as Float_t.
≅(a, b) = all(isapprox.(a, b; rtol=1e-5, atol=1e-7))

@testset "experiment styles and labels: validation and protocol (no ROOT)" begin
    # Unset: no scene, no style row, no label (the request and the plot are unchanged).
    @test isempty(scene_rows(GraphSpec()))
    @test ROOTPlots.experiment_label(GraphSpec()) === nothing
    @test isempty(scene_rows(Graph2DSpec()))

    # One switch: style=:CMS selects the style and the default CMS label.
    rows = scene_rows(GraphSpec(style=:CMS))
    @test rows[1] == "SCENE\t1" && rows[2] == "ROOTSTYLE\tCMS" && rows[end] == "ENDSCENE"
    @test label_row(GraphSpec(style=:CMS)) == ["EXPLABEL", "CMS", "inside_left", "Preliminary", "-", "0"]
    @test label_row(GraphSpec(style=:cms)) !== nothing                    # names are case-insensitive
    @test label_row(GraphSpec(style=:CMS, cms_label=false)) === nothing
    @test label_row(GraphSpec(style=:ATLAS)) === nothing                  # the ATLAS label is opt-in
    @test label_row(GraphSpec(style=:Plain)) === nothing

    l = ROOTPlots.experiment_label(GraphSpec(cms_label=(extra="Simulation", lumi="138 fb^{-1}",
                                                       energy="13 TeV", position=:outside)))
    @test (l.kind, l.extra, l.right, l.position) == ("CMS", "Simulation", "138 fb^{-1} (13 TeV)", :outside)
    @test ROOTPlots.experiment_label(GraphSpec(cms_label=(energy="13.6 TeV",))).right == "(13.6 TeV)"
    @test ROOTPlots.experiment_label(GraphSpec(cms_label=(lumi="1 fb^{-1}", extra=""))).extra == ""
    a = ROOTPlots.experiment_label(GraphSpec(atlas_label=(energy="13 TeV", lumi="140 fb^{-1}")))
    @test (a.kind, a.extra, a.lines, a.right) == ("ATLAS", "Internal", ["#sqrt{s} = 13 TeV, 140 fb^{-1}"], "")
    @test ROOTPlots.experiment_label(GraphSpec(atlas_label=true)).lines == String[]
    # An ATLAS label under style=:CMS replaces the default CMS label.
    @test ROOTPlots.experiment_label(GraphSpec(style=:CMS, atlas_label=true)).kind == "ATLAS"
    # 2-D/3-D views put the label above the frame unless a position is given.
    @test ROOTPlots.experiment_label(Graph2DSpec(style=:CMS)).position == :outside
    @test ROOTPlots.experiment_label(Graph3DSpec(cms_label=(position=:inside_right,))).position == :inside_right
    @test ROOTPlots.nd_label_fields(Graph2DSpec(style=:CMS)).cms_label == (position=:outside,)

    # Comparison panels: the label carries the panel's luminosity text and annotation.
    p = PanelSpec("x"; lumi="138 fb^{-1} (13 TeV)", annotation="ee channel", style=:CMS)
    l = ROOTPlots.experiment_label(p; lumi=p.lumi, annotation=p.annotation)
    @test (l.right, l.lines) == ("138 fb^{-1} (13 TeV)", ["ee channel"])
    p = PanelSpec("x"; lumi="140 fb^{-1}", atlas_label=true)
    @test ROOTPlots.experiment_label(p; lumi=p.lumi).lines == ["140 fb^{-1}"]
    legacy = PanelSpec("x", "Events", "", "", false, false, false, nothing, nothing, nothing,
                       AxisStyle(), AxisStyle(), AxisStyle(), CanvasStyle())
    @test legacy.style === nothing && legacy.cms_label === nothing

    for bad in (GraphSpec(style=Symbol("no space")), GraphSpec(style=Symbol("")),
                GraphSpec(cms_label=(colour="red",)), GraphSpec(cms_label=(position=:top,)),
                GraphSpec(cms_label=(extra=1,)), GraphSpec(cms_label=true, atlas_label=true),
                GraphSpec(atlas_label=(extra="Internal",)), GraphSpec(cms_label=(extra="a\tb",)))
        @test_throws ArgumentError ROOTPlots.validate_style_and_label(bad)
    end
    @test_throws ArgumentError ROOTPlots.validate_nd_spec(Graph2DSpec(cms_label=(foo=1,)))
    @test_throws ArgumentError PlotFigure([PlotPanel(Hist1D([0., 1.], [1.], [1.]), GraphSpec())]; style=Symbol("a b"))

    # gROOT->SetStyle is process-global: panels and figures must agree (checked before ROOT is needed).
    h = Hist1D([0., 1., 2.], [1., 2.], [1., 2.])
    mixed = PlotFigure([PlotPanel(h, GraphSpec(style=:CMS)), PlotPanel(h, GraphSpec(style=:ATLAS))])
    @test_throws ArgumentError render_root(mixed, joinpath(tempdir(), "never-written"))
    @test_throws ArgumentError render_root(PlotFigure([PlotPanel(h, GraphSpec(style=:ATLAS))]; style=:CMS),
                                           joinpath(tempdir(), "never-written"))
    # The figure style reaches every panel spec without a style.
    styled = ROOTPlots.styled_figure(PlotFigure([PlotPanel(h, GraphSpec()), PlotFigure([PlotPanel(h, Graph2DSpec())])]; style=:CMS))
    @test all(p -> p.spec.style == :CMS, ROOTPlots.figure_leaves(styled))
    cmp = comparison(h, [Component("MC", h, "#336699")])
    @test_throws ArgumentError render_root([cmp, cmp], [PanelSpec("x"; style=:CMS), PanelSpec("x"; style=:Plain)],
                                           joinpath(tempdir(), "never-written"))
end

if TEST_NATIVE_ROOT
@testset "experiment styles reach ROOT (native)" begin
    @test issubset(["ATLAS", "CMS", "Modern", "Plain"], root_styles())
    s = GraphSeries(collect(1.0:10.0), collect(11.0:20.0); yerrors=fill(1.0, 10), label="Data")

    # Unset: ROOT's default style, no label.
    c = only(with_axis_dump(dir -> render_root([s], GraphSpec(), joinpath(dir, "plain"))))
    @test c.style[].name == "Modern"
    @test isempty(only(framed(c)).latex)

    # One switch: the CMS TStyle values and the CMS label in the pad.
    spec = GraphSpec(xlabel="m [GeV]", ylabel="Events", style=:CMS, cms_label=(lumi="138 fb^{-1}", energy="13 TeV"))
    c = only(with_axis_dump(dir -> render_root([s], spec, joinpath(dir, "cms"))))
    st = c.style[]
    @test st.name == "CMS" && st.margins ≅ (0.16, 0.04, 0.13, 0.08)
    @test (st.label_font, st.title_font) == (42, 42) && (st.label_size, st.title_size) ≅ (0.05, 0.06)
    @test (st.tickx, st.ticky, st.optstat) == (1, 1, 0)
    pad = only(framed(c))
    @test pad.margins ≅ (0.16, 0.04, 0.13, 0.08)
    @test (pad.tickx, pad.ticky) == (1, 1)
    @test pad.axes["x"].label_size ≅ 0.05 && pad.axes["y"].title_size ≅ 0.06 && pad.axes["x"].label_font == 42
    cms = only(latex(pad, "CMS")); pre = only(latex(pad, "Preliminary")); lumi = only(latex(pad, "138 fb^{-1} (13 TeV)"))
    @test (cms.font, cms.align, cms.ndc) == (61, 13, true) && cms.size ≅ 0.75 * 0.08
    @test cms.x ≅ 0.16 + 0.045 * 0.80 && cms.y ≅ 1 - 0.08 - 0.035 * 0.79
    @test (pre.font, pre.align) == (52, 13) && pre.size ≅ 0.76 * 0.75 * 0.08 && pre.y < cms.y
    @test (lumi.font, lumi.align) == (42, 31) && lumi.size ≅ 0.6 * 0.08 && lumi.x ≅ 0.96 && lumi.y > 0.92

    # Explicit options still win over the style.
    spec = GraphSpec(style=:CMS, canvas=CanvasStyle(left_margin=0.2, ticks_y=false), xaxis=AxisStyle(label_size=0.03))
    pad = only(framed(only(with_axis_dump(dir -> render_root([s], spec, joinpath(dir, "override"))))))
    @test pad.margins[1] ≅ 0.2 && pad.ticky == 0 && pad.axes["x"].label_size ≅ 0.03 && pad.axes["y"].label_size ≅ 0.05
    @test only(latex(pad, "CMS")).x ≅ 0.2 + 0.045 * 0.76

    # Outside position: "CMS" and the extra text on one line above the frame.
    spec = GraphSpec(style=:CMS, cms_label=(extra="Simulation", position=:outside))
    pad = only(framed(only(with_axis_dump(dir -> render_root([s], spec, joinpath(dir, "outside"))))))
    cms = only(latex(pad, "CMS")); sim = only(latex(pad, "Simulation"))
    @test cms.align == 11 && sim.align == 11 && cms.y ≅ sim.y && cms.y > 0.92 && sim.x > cms.x && sim.font == 52

    # ROOT's built-in ATLAS style and the ATLAS label.
    spec = GraphSpec(style=:ATLAS, atlas_label=(text="Internal", energy="13 TeV", lumi="140 fb^{-1}"))
    c = only(with_axis_dump(dir -> render_root([s], spec, joinpath(dir, "atlas"))))
    @test c.style[].name == "ATLAS"
    pad = only(framed(c))
    atlas = only(latex(pad, "ATLAS")); internal = only(latex(pad, "Internal"))
    line = only(latex(pad, "#sqrt{s} = 13 TeV, 140 fb^{-1}"))
    @test (atlas.font, internal.font, line.font) == (72, 42, 42)
    @test atlas.y ≅ internal.y && internal.x > atlas.x && line.y < atlas.y && line.x ≅ atlas.x

    # Unknown names fail in the native backend with the list of available styles.
    exe = ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
    mktempdir() do dir
        bad = GraphSpec(style=:NoSuchStyle)
        input = joinpath(dir, "graph.tsv")
        ROOTPlots.write_graph_spec(input, [s], bad, ROOTPlots.graph_limits([s], bad))
        err = IOBuffer()
        proc = run(pipeline(ignorestatus(Cmd([exe, input, joinpath(dir, "out")])); stderr=err))
        message = String(take!(err))
        @test proc.exitcode != 0
        @test occursin("unknown ROOT style \"NoSuchStyle\"", message) && occursin("ATLAS", message) &&
              occursin("CMS", message) && occursin("Plain", message)
    end

    # Hist1D + fit, fit/pull view: the label sits in the upper pad, sized as in one pad.
    # Sum of four uniform grids: a nearly Gaussian, deterministic sample.
    values = vec([a + b + c + d for a in -1:0.25:1, b in -1:0.25:1, c in -1:0.25:1, d in -1:0.25:1]) .+ 0.1
    h = histogram(values, collect(range(-4.0, 4.0; length=17))).hist
    fit = fit_root(h, :gaus; initial=[1500.0, 0.0, 1.2])
    @test fit.status == 0 && fit.valid
    c = only(with_axis_dump(dir -> render_root(h, fit, GraphSpec(style=:CMS), joinpath(dir, "fitpull"))))
    top = only(p for p in c.pads if p.name == "fit_top"); bottom = only(p for p in c.pads if p.name == "fit_bottom")
    @test top.margins[4] ≅ 0.08 / 0.70 && bottom.margins[3] ≅ 0.13 / 0.30 && top.margins[1] ≅ 0.16
    @test only(latex(top, "CMS")).size ≅ 0.75 * 0.08 * 600 / 420 && isempty(latex(bottom, "CMS"))
    @test top.axes["y"].label_size ≅ 0.05 * 600 / 420 && bottom.axes["x"].label_size ≅ 0.05 * 600 / 180

    # Comparison with ratio panel: luminosity text and annotation go through the label.
    cmp = comparison(h, [Component("MC", h, "#336699")])
    spec = PanelSpec("x"; lumi="138 fb^{-1} (13 TeV)", annotation="ee", style=:CMS)
    c = only(with_axis_dump(dir -> render_root([cmp], [spec], joinpath(dir, "ratio"))))
    @test c.style[].name == "CMS"
    top = only(p for p in c.pads if p.name == "top0")
    @test top.margins[4] ≅ 0.08 / 0.70 && top.margins[1] ≅ 0.16
    @test length(latex(top, "CMS")) == 1 && length(latex(top, "138 fb^{-1} (13 TeV)")) == 1 && length(latex(top, "ee")) == 1

    # Hist2D COLZ: style margins, palette room on the right, label above the frame.
    xy = histogram((collect(range(-1.0, 1.0; length=200)), collect(range(-1.0, 1.0; length=200)) .^ 2),
                   (collect(range(-1.0, 1.0; length=11)), collect(range(0.0, 1.0; length=11)))).hist
    c = only(with_axis_dump(dir -> render_root(xy, Graph2DSpec(style=:CMS), joinpath(dir, "colz"); draw="COLZ")))
    pad = only(framed(c))
    @test pad.margins[1] ≅ 0.16 && pad.margins[2] >= 0.18 && pad.margins[3] ≅ 0.13
    @test only(latex(pad, "CMS")).align == 11 && only(latex(pad, "CMS")).y > 1 - pad.margins[4]

    # 2x2 figure: one style for the canvas, a label in each pad that asks for one.
    figure = PlotFigure([PlotPanel(h, GraphSpec()), PlotPanel(xy, Graph2DSpec()),
                         PlotPanel(h, GraphSpec(cms_label=(extra="Simulation",))),
                         PlotPanel([s], GraphSpec(cms_label=false))]; layout=(2, 2), style=:CMS)
    c = only(with_axis_dump(dir -> render_root(figure, joinpath(dir, "figure"))))
    pads = framed(c)
    @test c.style[].name == "CMS" && length(pads) == 4
    @test [length(latex(p, "CMS")) for p in pads] == [1, 1, 1, 0]
    @test length(latex(pads[3], "Simulation")) == 1
    @test only(latex(pads[1], "CMS")).size ≅ 0.75 * pads[1].margins[4]

    # A titled COLZ palette fits its pad under every style: ROOT's TPaletteAxis
    # sits 0.005-0.05 pad widths right of the frame, its labels, "x10^n" exponent
    # and z title further right. The dump gives the palette and z-title extents
    # (pad NDC); the PNG gives the dark pixels that ROOT actually painted.
    @testset "titled palette fits: $(something(style, :default))" for style in (nothing, :CMS, :ATLAS)
        w = [round(236 * exp(-((i - 5.5)^2 + (j - 4.5)^2) / 8)) for i in 1:10, j in 1:8]
        small = Hist2D(collect(0.0:10.0), collect(0.0:8.0), w, w)
        large = Hist2D(collect(0.0:10.0), collect(0.0:8.0), 1e4 .* w, 1e8 .* w)   # labels x10^{3}
        label = style == :CMS ? (cms_label=(lumi="138 fb^{-1}", energy="13 TeV"),) :
                style == :ATLAS ? (atlas_label=(energy="13 TeV", lumi="140 fb^{-1}"),) : (;)
        right_text = style == :CMS ? "138 fb^{-1} (13 TeV)" : "#sqrt{s} = 13 TeV, 140 fb^{-1}"
        # Descenders and parentheses reach furthest right in a rotated title.
        ztitle = "Weighted yield (pb)"
        spec = Graph2DSpec(; xlabel="x", ylabel="y", zlabel=ztitle, style, label...)
        fits(pad) = (p = pad.palette[];
                     p !== nothing && p.x1 > 1 - pad.margins[2] && p.title_left > p.x2 + pad.axes["z"].label_offset &&
                     p.title_right <= 0.99)
        mktempdir() do keep
            for (name, h) in ("small" => small, "large" => large)
                png = joinpath(keep, name * ".png")
                c = only(with_axis_dump(dir -> cp(render_root(h, spec, joinpath(dir, name); draw="COLZ").png, png)))
                pad = only(framed(c)); p = pad.palette[]
                @test fits(pad)
                # ROOT's default style keeps its look while the labels are narrow.
                style === nothing && name == "small" && @test pad.margins[2] ≅ 0.2 && pad.axes["z"].title_offset ≅ 1.3
                if style !== nothing
                    # The experiment label stays above the frame, left of the palette and its exponent.
                    @test only(latex(pad, right_text)).x <= p.x1 && only(latex(pad, right_text)).align == 31
                    @test only(latex(pad, string(style))).y > p.y2
                end
                # Pixels: nothing touches the right canvas edge (a clipped title would), and
                # nothing above the palette touches the top edge (a clipped exponent would).
                ink = ROOTPlots.image_ink_bounds(png)
                @test 0 < ink.xmax <= ink.width - 9
                top = ROOTPlots.image_ink_bounds(png; region=(round(Int, p.x1 * ink.width), 0, ink.width - 1, ink.height - 1))
                @test top.ymin >= 2
                # A small style top margin (ATLAS: 0.05) grows for the exponent above the
                # palette and for the outside label (both here).
                style == :ATLAS && @test pad.margins[4] > 0.05 + 1e-6
            end
            # Figure pads: the left pad ends inside the canvas, the right one at its edge.
            fig = PlotFigure([PlotPanel(large, Graph2DSpec(; xlabel="x", ylabel="y", zlabel=ztitle, label...); draw="COLZ"),
                              PlotPanel(small, Graph2DSpec(; xlabel="x", ylabel="y", zlabel=ztitle, label...); draw="COLZ")];
                             layout=(1, 2), size=(1200, 450), style)
            png = joinpath(keep, "figure.png")
            c = only(with_axis_dump(dir -> cp(render_root(fig, joinpath(dir, "figure")).png, png)))
            pads = framed(c)
            @test length(pads) == 2 && all(fits, pads)
            full = ROOTPlots.image_ink_bounds(png)
            W, H = full.width, full.height
            for pad in pads
                x0, y0 = round(Int, pad.place[1] * W), round(Int, (1 - pad.place[2] - pad.place[4]) * H)
                x1, y1 = round(Int, (pad.place[1] + pad.place[3]) * W) - 1, round(Int, (1 - pad.place[2]) * H) - 1
                ink = ROOTPlots.image_ink_bounds(png; region=(x0, y0, x1, y1))
                @test 0 < ink.xmax <= x1 - 8
                top = ROOTPlots.image_ink_bounds(png; region=(x0 + round(Int, pad.palette[].x1 * (x1 - x0 + 1)), y0, x1, y1))
                @test top.ymin >= y0 + 2
            end
        end
    end

    # Explicit choices win: a fixed right margin and z title offset are kept;
    # a larger explicit z title size still gets room.
    w = [Float64(i * j) for i in 1:10, j in 1:8]
    h = Hist2D(collect(0.0:10.0), collect(0.0:8.0), w, w)
    spec = Graph2DSpec(zlabel="Events", style=:CMS, canvas=CanvasStyle(right_margin=0.3), zaxis=AxisStyle(title_offset=2.0))
    pad = only(framed(only(with_axis_dump(dir -> render_root(h, spec, joinpath(dir, "fixed"); draw="COLZ")))))
    @test pad.margins[2] ≅ 0.3 && pad.axes["z"].title_offset ≅ 2.0
    spec = Graph2DSpec(zlabel="Events", style=:CMS, zaxis=AxisStyle(title_size=0.08))
    pad = only(framed(only(with_axis_dump(dir -> render_root(h, spec, joinpath(dir, "big"); draw="COLZ")))))
    @test pad.axes["z"].title_size ≅ 0.08 && pad.palette[].title_right <= 0.99 && pad.margins[2] > 0.2
end

    # Label texts above the frame end inside their pad with white pixels at its
    # top edge, for every view, under small style top margins (ATLAS: 0.05).
    # Above-frame texts are baseline-aligned (vertical alignment 1): their TLatex
    # box spans y to y + height (pad NDC, from the dump).
    @testset "outside labels fit inside the pad: $kind" for kind in (:ATLAS, :CMS)
        w = [round(236 * exp(-((i - 5.5)^2 + (j - 4.5)^2) / 8)) for i in 1:10, j in 1:8]
        small = Hist2D(collect(0.0:10.0), collect(0.0:8.0), w, w)
        large = Hist2D(collect(0.0:10.0), collect(0.0:8.0), 1e4 .* w, 1e8 .* w)   # palette x10^{3}
        outside(position=:outside) = kind == :ATLAS ?
            (atlas_label=(energy="13 TeV", lumi="140 fb^{-1}", position=position),) :
            (cms_label=(lumi="138 fb^{-1}", energy="13 TeV", position=position),)
        first_text, line_text = kind == :ATLAS ? ("ATLAS", "#sqrt{s} = 13 TeV, 140 fb^{-1}") : ("CMS", "138 fb^{-1} (13 TeV)")
        above(pad) = filter(t -> t.align % 10 == 1, pad.latex)
        # Dump: every text above the frame starts above it (frame=true) and its box
        # ends inside the pad. Pixels: at least 3 white rows at the pad's top edge.
        function label_fits(pad, png; nabove=1, frame=true)
            texts = above(pad)
            length(texts) >= nabove || return false
            all(t -> (!frame || t.y >= 1 - pad.margins[4] - 1e-9) && t.y + t.height <= 1, texts) || return false
            full = ROOTPlots.image_ink_bounds(png)
            x0, y0 = round(Int, pad.place[1] * full.width), round(Int, (1 - pad.place[2] - pad.place[4]) * full.height)
            x1, y1 = round(Int, (pad.place[1] + pad.place[3]) * full.width) - 1, round(Int, (1 - pad.place[2]) * full.height) - 1
            return ROOTPlots.image_ink_bounds(png; region=(x0, y0, x1, y1)).ymin >= y0 + 3
        end
        # Render, keep the PNG while its pads are checked.
        fits_in(f, name, select; nabove=1, frame=true) = mktempdir() do keep
            png = joinpath(keep, name * ".png")
            c = only(with_axis_dump(dir -> cp(f(joinpath(dir, name)).png, png)))
            pads = select(c)
            (c, pads, [label_fits(p, png; nabove, frame) for p in pads])
        end

        # 2-D COLZ, default (outside) label, with and without the palette exponent.
        for (name, h) in ("small" => small, "large" => large)
            c, pads, ok = fits_in(stem -> render_root(h, Graph2DSpec(; xlabel="x", ylabel="y", zlabel="Events",
                                                                      style=kind, outside()...), stem; draw="COLZ"),
                                  "colz_$name", c -> framed(c); nabove=2)
            pad = only(pads)
            @test only(ok)
            # One row above the frame, left of the palette (and its exponent).
            head = only(latex(pad, first_text)); line = only(latex(pad, line_text))
            @test head.y ≅ line.y && line.align == 31 && line.x <= pad.palette[].x1
            # ATLAS: the 0.05 style top margin grows for the label.
            kind == :ATLAS && @test pad.margins[4] > 0.05 + 1e-3
            kind == :CMS && @test pad.margins[4] ≅ 0.08
        end

        # Figure pads: 1x2 (one row per pad) and 1x3 (too narrow: ATLAS stacks the lines).
        for n in (2, 3)
            fig = PlotFigure([PlotPanel(n == 2 && i == 1 ? large : small,
                                        Graph2DSpec(; xlabel="x", ylabel="y", zlabel="Events", outside()...); draw="COLZ")
                              for i in 1:n]; layout=(1, n), size=(1200, 450), style=kind)
            c, pads, ok = fits_in(stem -> render_root(fig, stem), "figure_$n", c -> framed(c); nabove=2)
            @test length(pads) == n && all(ok)
            for pad in pads
                head = only(latex(pad, first_text)); line = only(latex(pad, line_text))
                if kind == :ATLAS && n == 3
                    @test line.align == 11 && line.x ≅ head.x && line.y < head.y
                else
                    @test head.y ≅ line.y && line.x <= pad.palette[].x1
                end
            end
        end

        # 1-D views: outside and inside labels on one canvas, the fit/pull and
        # comparison upper pads (outside), and a figure pad.
        s = GraphSeries(collect(1.0:10.0), collect(11.0:20.0); yerrors=fill(1.0, 10), label="Data")
        for position in (:outside, :inside_left, :inside_right)
            c, pads, ok = fits_in(stem -> render_root([s], GraphSpec(; style=kind, outside(position)...), stem),
                                  "graph_$position", c -> framed(c);
                                  nabove=position == :outside ? 2 : kind == :CMS ? 1 : 0)
            @test only(ok)
        end
        values = vec([a + b + c + d for a in -1:0.25:1, b in -1:0.25:1, c in -1:0.25:1, d in -1:0.25:1]) .+ 0.1
        h = histogram(values, collect(range(-4.0, 4.0; length=17))).hist
        fit = fit_root(h, :gaus; initial=[1500.0, 0.0, 1.2])
        c, pads, ok = fits_in(stem -> render_root(h, fit, GraphSpec(; style=kind, outside()...), stem),
                              "fitpull", c -> [only(p for p in c.pads if p.name == "fit_top")]; nabove=2)
        @test only(ok)
        cmp = comparison(h, [Component("MC", h, "#336699")])
        c, pads, ok = fits_in(stem -> render_root([cmp], [PanelSpec("x"; style=kind, outside()...)], stem),
                              "ratio", c -> [only(p for p in c.pads if p.name == "top0")]; nabove=2)
        @test only(ok)
        fig = PlotFigure([PlotPanel(h, GraphSpec(; outside()...)), PlotPanel(small, Graph2DSpec(; outside()...); draw="COLZ")];
                         layout=(1, 2), size=(1000, 400), style=kind)
        c, pads, ok = fits_in(stem -> render_root(fig, stem), "figure_1d", c -> framed(c); nabove=2)
        @test length(pads) == 2 && all(ok)

        # An explicit top margin is kept; the texts move down to stay inside the pad.
        spec = Graph2DSpec(; style=kind, canvas=CanvasStyle(top_margin=0.03), outside()...)
        c, pads, ok = fits_in(stem -> render_root(small, spec, stem; draw="COLZ"), "fixed", c -> framed(c);
                              nabove=2, frame=false)
        @test only(pads).margins[4] ≅ 0.03 && only(ok)
    end
end

# Two-pad views (comparison ratio pad, fit/pull): text has one absolute (pixel)
# size in both pads, the two y titles line up, and labels at the frame edges
# (1.5 and 0.5 on a 0.5-1.5 ratio axis, 0 on a linear upper axis) end inside
# their pads. Pixel sizes come from the dump (ROOT's sizes are fractions of the
# smaller pad dimension, or pixels for precision-3 fonts) and the PNG size.
if TEST_NATIVE_ROOT
"""Pixel geometry of the two pads of a dumped two-pad canvas rendered to png."""
function two_pad_geometry(c, png, upper, lower)
    full = ROOTPlots.image_ink_bounds(png)
    W, H = full.width, full.height
    function geometry(p)
        w, h = p.place[3] * W, p.place[4] * H
        px(size, font) = font % 10 == 3 ? size : size * min(w, h)
        unit(a, length) = a.title_font % 10 == 3 ? a.title_size : a.title_size * length
        x0, y0 = round(Int, p.place[1] * W), round(Int, (1 - p.place[2] - p.place[4]) * H)
        x1, y1 = round(Int, (p.place[1] + p.place[3]) * W) - 1, round(Int, (1 - p.place[2]) * H) - 1
        (; pad=p, w, h, x0, y0, x1, y1, frame_left=x0 + round(Int, p.margins[1] * w),
         frame_top=y0 + round(Int, p.margins[4] * h), frame_bottom=y1 - round(Int, p.margins[3] * h),
         y_label=px(p.axes["y"].label_size, p.axes["y"].label_font), y_title=px(p.axes["y"].title_size, p.axes["y"].title_font),
         x_label=px(p.axes["x"].label_size, p.axes["x"].label_font), x_title=px(p.axes["x"].title_size, p.axes["x"].title_font),
         # TGaxis: a y title 1.6*size*offset pad widths left of the axis
         y_title_distance=1.6 * p.axes["y"].title_offset * unit(p.axes["y"], w))
    end
    return geometry(only(q for q in c.pads if q.name == upper)), geometry(only(q for q in c.pads if q.name == lower))
end
"""Render with f(stem), then run check(c, upper, lower, png) while the PNG exists."""
two_pad_render(check, f, name, upper, lower) = mktempdir() do keep
    png = joinpath(keep, name * ".png")
    c = only(with_axis_dump(dir -> cp(f(joinpath(dir, name)).png, png)))
    check(c, two_pad_geometry(c, png, upper, lower)..., png)
end
# Equal pixel sizes and aligned titles; room for half a y label at both inner frame edges;
# ink: the edge labels (and anything else) stay at least 2 pixels inside the pads.
function two_pad_ok(t, b, png)
    ok = isapprox(b.y_label, t.y_label; rtol=0.01) && isapprox(b.x_label, t.y_label; rtol=0.01) &&
         isapprox(b.x_title, t.y_title; rtol=0.01) && isapprox(b.y_title, t.y_title; rtol=0.01) &&
         abs(b.y_title_distance - t.y_title_distance) <= 1
    ok &= b.pad.margins[4] * b.h >= 0.5 * b.y_label + 2 && t.pad.margins[3] * t.h >= 0.5 * t.y_label + 2
    # Left of the lower frame, from the pad top to below the frame top: the top label (1.5) is inside.
    top_label = ROOTPlots.image_ink_bounds(png; region=(b.x0, b.y0, b.frame_left - 2, b.frame_top + round(Int, b.y_label)))
    ok &= top_label.ymin >= b.y0 + 2
    # Left of the upper frame, its lowest labels end above the pad's bottom edge.
    low_label = ROOTPlots.image_ink_bounds(png; region=(t.x0, t.frame_bottom - round(Int, t.y_label), t.frame_left - 2, t.y1))
    ok &= low_label.ymax == -1 || low_label.ymax <= t.y1 - 2
    # Nothing in the upper pad (exponent, labels, titles) reaches its top edge.
    ok &= ROOTPlots.image_ink_bounds(png; region=(t.x0, t.y0, t.x1, t.y1)).ymin >= t.y0 + 2
    return ok
end

@testset "two-pad text and edge labels" begin
    # A CMS-like data/MC comparison: falling stacked prediction, counts in the thousands.
    edges = collect(35.0:10.0:205.0)
    zj = [90, 55, 20, 13, 12, 6.5, 5.5, 3.5, 2.5, 1.2, 2.0, 1.0, 0.6, 0.7, 0.6, 0.6, 2.5]
    total = [2800, 2200, 1800, 1400, 1100, 900, 700, 520, 380, 290, 220, 170, 140, 110, 85, 60, 260]
    mk(v) = Hist1D(edges, Float64.(v), 0.02 .* v)
    cmp = comparison(Hist1D(edges, round.(0.9 .* total), round.(0.9 .* total)),
                     [Component("Z + jets", mk(zj), "#1f77b4"), Component("t#bar{t}", mk(total .- zj), "#c00020")];
                     uncertainty_label="MC stat. unc.")
    values = vec([a + b + c + d for a in -1:0.25:1, b in -1:0.25:1, c in -1:0.25:1, d in -1:0.25:1]) .+ 0.1
    h = histogram(values, collect(range(-4.0, 4.0; length=17))).hist
    fit = fit_root(h, :gaus; initial=[1500.0, 0.0, 1.2])
    for style in (nothing, :CMS, :ATLAS), logy in (false, true)
        label = style == :CMS ? (cms_label=(extra="tutorial", lumi="1.21 fb^{-1}", energy="13.6 TeV"),) :
                style == :ATLAS ? (atlas_label=(energy="13.6 TeV", lumi="1.21 fb^{-1}"),) : (;)
        spec = PanelSpec("Lepton p_{T} (GeV)"; ylabel="Events / 10 GeV", logy, ratio_limits=(0.5, 1.5), style, label...)
        two_pad_render(stem -> render_root([cmp], [spec], stem), "ratio", "top0", "bottom0") do c, t, b, png
            @test two_pad_ok(t, b, png)
            @test b.pad.axes["y"].title == "Data / pred."
        end
        two_pad_render(stem -> render_root(h, fit, GraphSpec(; xlabel="x", ylabel="Entries", logy, style), stem),
                       "fitpull", "fit_top", "fit_bottom") do c, t, b, png
            @test two_pad_ok(t, b, png)
        end
    end
    # The default style's title sizes are kept (upper pad as before; the ratio pad follows it).
    two_pad_render(stem -> render_root([cmp], [PanelSpec("x")], stem), "plain", "top0", "bottom0") do c, t, b, png
        @test t.pad.axes["y"].title_size ≅ 0.035 && t.pad.axes["y"].label_size ≅ 0.038 && t.pad.axes["y"].title_offset ≅ 1.55
    end
    # Under CMS the upper y title keeps the style's one-pad size (0.06 of the canvas width).
    two_pad_render(stem -> render_root([cmp], [PanelSpec("x"; style=:CMS)], stem), "cms", "top0", "bottom0") do c, t, b, png
        @test isapprox(t.y_title, 0.06 * t.w; rtol=0.01) && isapprox(t.y_label, 0.05 * t.w; rtol=0.01)
    end

    # Explicit sizes win, in ROOT's units for their own pad; unset lower-pad sizes follow
    # the upper y axis (a precision-3 font gives pixels).
    spec = PanelSpec("x"; style=:CMS, yaxis=AxisStyle(title_size=0.035, label_size=0.035))
    two_pad_render(stem -> render_root([cmp], [spec], stem), "small", "top0", "bottom0") do c, t, b, png
        @test t.pad.axes["y"].title_size ≅ 0.035 && two_pad_ok(t, b, png)
        @test isapprox(b.y_label, 0.035 * min(t.w, t.h); rtol=0.01)
    end
    spec = PanelSpec("x"; yaxis=AxisStyle(title_font=43, title_size=24), ratio_yaxis=AxisStyle(label_size=0.12))
    two_pad_render(stem -> render_root([cmp], [spec], stem), "explicit", "top0", "bottom0") do c, t, b, png
        @test t.pad.axes["y"].title_size ≅ 24 && b.pad.axes["y"].label_size ≅ 0.12 && b.pad.axes["y"].title_font == 42
        @test isapprox(b.y_title, 24; rtol=0.01) && isapprox(b.x_title, 24; rtol=0.01)
        @test abs(b.y_title_distance - t.y_title_distance) <= 1
    end
    spec = PanelSpec("x"; ratio_yaxis=AxisStyle(title_offset=2.0))
    two_pad_render(stem -> render_root([cmp], [spec], stem), "offset", "top0", "bottom0") do c, t, b, png
        @test b.pad.axes["y"].title_offset ≅ 2.0
    end
    # ROOT's automatic y title offset (the default style's 0, here with a label only) would
    # differ between the pads; the titles still line up.
    two_pad_render(stem -> render_root([cmp], [PanelSpec("x"; cms_label=true)], stem), "auto", "top0", "bottom0") do c, t, b, png
        @test t.pad.axes["y"].title_offset > 0 && two_pad_ok(t, b, png)
    end
end
end
