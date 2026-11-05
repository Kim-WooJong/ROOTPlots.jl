# Unicode labels -> ROOT TLatex (root_latex) and its uniform use for every text sent to ROOT.
using ROOTPlots: root_latex, root_text, check_ascii_request

@testset "Unicode labels: TLatex translation (no ROOT)" begin
    # ASCII, including ROOT's own #-syntax, is passed through exactly.
    for s in ("", "x", "p_{T} [GeV]", "#alpha#beta", "#sqrt{s} = 13 TeV", "a{}b ^{} _x", "#font[12]{l}")
        @test root_latex(s) === s
    end
    @test root_latex("σ_m [GeV/c²] → Δφ") == "#sigma_m [GeV/c^{2}] #rightarrow #Delta#phi"
    @test root_latex("α β") == "#alpha #beta"
    lower = "αβγδεζηθικλμνξοπρστυφχψω"
    upper = "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ"
    names = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
             "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi",
             "chi", "psi", "omega"]
    @test root_latex(lower) == join("#" .* names)
    @test root_latex(upper) == join("#" .* uppercasefirst.(names))
    @test root_latex("ϵϑϖςϕϒµ") == "#varepsilon#vartheta#varomega#varsigma#varphi#varUpsilon#mu"
    @test root_latex("Ω") == "#Omega"                       # U+2126 ohm sign (NFD -> Omega)
    # Superscript/subscript runs become one group.
    @test root_latex("x⁻¹ cm⁻² s⁻¹") == "x^{-1} cm^{-2} s^{-1}"
    @test root_latex("H₂O, x₁₂, (n)ⁱ") == "H_{2}O, x_{12}, (n)^{i}"
    @test root_latex("10⁺³") == "10^{+3}"
    @test root_latex("π²") == "#pi^{2}"
    # Symbols, arrows and typography.
    @test root_latex("± ∓ × ÷ · ° − ≤ ≥ ≠ ≈ ≡ ∝ ∞ ∂ ∇ √ ∑ ∫") ==
          "#pm #mp #times #divide #upoint #circ #minus #leq #geq #neq #approx #equiv #propto #infty #partial #nabla #surd #sum #int"
    @test root_latex("← → ↑ ↓ ↔ ⇐ ⇒ ⇔") ==
          "#leftarrow #rightarrow #uparrow #downarrow #leftrightarrow #Leftarrow #Rightarrow #Leftrightarrow"
    @test root_latex("ℓ⁺ℓ⁻ ℏ ħ Å Å å") == "#font[12]{l}^{+}#font[12]{l}^{-} #hbar #hbar #AA #AA #aa"
    @test root_latex("1–2 “a” it’s") == "1-2 \"a\" it's"
    @test root_latex("5\u00a0°C") == "5 #circC"
    # Accents (precomposed or combining) use TLatex accents.
    @test root_latex("é ü ñ x̄ ŷ ẋ p⃗ θ̂") == "#acute{e} #ddot{u} #tilde{n} #bar{x} #hat{y} #dot{x} #vec{p} #hat{#theta}"
    # A command that a following letter would extend is ended by ROOT's "^{}".
    @test root_latex("∈t ⊂eq ∈ x") == "#in^{}t #subset^{}eq #in x"
    @test root_latex("αx μm") == "#alphax #mum"      # TLatex reads Greek names by prefix
    # Unsupported characters fail clearly (no mojibake), naming the character.
    for (bad, char) in (("\uc9c8\ub7c9 [GeV]", "\uc9c8"), ("ε ☃", "☃"), ("x̊", "̊"), ("\u65e5\u672c", "\u65e5"))
        err = try root_latex(bad); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("ROOT TLatex cannot draw '$char'", err.msg) && occursin("U+", err.msg)
        @test occursin("#-syntax", err.msg)
    end
    @test_throws ArgumentError root_text("a\tb")
    @test_throws ArgumentError root_text("a\0b")
    @test all(isascii(root_latex(string(c))) for c in keys(ROOTPlots.ROOT_LATEX_SYMBOLS))
end

# The caller's ASCII text (ROOT #-syntax included) is copied verbatim; the "^{}"
# separator ends only commands produced by the translation itself.
const MIXED_ROOT_LABELS = [
    "#sigma_{m} [GeV/c²]" => "#sigma_{m} [GeV/c^{2}]",
    "Δ#phi" => "#Delta#phi",
    "#splitline{α}{β}" => "#splitline{#alpha}{#beta}",
    "#it{p}_{T} [μm]" => "#it{p}_{T} [#mum]",
    "x #in A, α" => "x #in A, #alpha",
    "#int α" => "#int #alpha",
    "#infty · α" => "#infty #upoint #alpha",
    "#subset ⊂ #supset" => "#subset #subset #supset",
    "#in∈t" => "#in#in^{}t",
    "α#" => "#alpha#",
    "# α" => "# #alpha",
    "##α" => "###alpha",
    "#bf{μ} #in [0, 1]" => "#bf{#mu} #in [0, 1]",
    "#chi^{2}/ndf ≈ 1" => "#chi^{2}/ndf #approx 1",
    "#sqrt{s} = 13 TeV, ∫L = 138 fb⁻¹" => "#sqrt{s} = 13 TeV, #intL = 138 fb^{-1}",
]

@testset "Unicode mixed with ROOT #-syntax (no ROOT)" begin
    for (label, expected) in MIXED_ROOT_LABELS
        @test root_latex(label) == expected
        # Every ASCII run of the label appears unchanged and in order.
        out, at = root_latex(label), 1
        for run in eachmatch(r"[\x00-\x7f]+", label)
            r = findnext(run.match, out, at)
            @test r !== nothing
            r === nothing || (at = last(r) + 1)
        end
    end
    # The prefix guard for translated #in / #subset / #supset is unchanged.
    @test root_latex("∈t ⊂eq ⊃eq ∈f") == "#in^{}t #subset^{}eq #supset^{}eq #in^{}f"
    @test root_latex("∈ t") == "#in t" && root_latex("∈²") == "#in^{2}" && root_latex("∈1") == "#in1"
    # needs_latex_separator: only complete command names that a letter extends.
    nls = ROOTPlots.needs_latex_separator
    @test nls("#in") && nls("#subset") && nls("#supset") && nls("#i") && nls("#inf")
    @test !any(nls, ("#", "", "in", "#1", "#in{", "#int", "#infty", "#subseteq", "#alpha", "#font[12]{l}", "^{2}", "#in^{}"))
    # Idempotence: the translation is ASCII, which is passed through unchanged.
    battery = [first.(MIXED_ROOT_LABELS); "σ_m [GeV/c²] → Δφ"; "∈t ⊂eq"; "é ü x̄ p⃗"; "ℓ⁺ℓ⁻ Å"; "x⁻¹ H₂O"]
    for label in battery
        once = root_latex(label)
        @test isascii(once) && root_latex(once) === once && root_text(once) == once
    end
end

@testset "GraphSeries labels: drawing translates, fitting ignores (no ROOT)" begin
    g = GraphSeries([1.0, 2.0], [1.0, 2.0]; yerrors=[0.1, 0.1], label="\uce21\uc815 \U1F600")
    @test g.label == "\uce21\uc815 \U1F600"                       # stored verbatim
    @test ROOTPlots.validate_graph(g; label=false) === g             # fit validation: numbers only
    @test_throws ArgumentError ROOTPlots.validate_graph(g)            # drawing: TLatex cannot draw it
    @test_throws ArgumentError ROOTPlots.graph_limits([g], GraphSpec())
    # Numerical checks still apply without the label.
    bad = GraphSeries(g.x, g.y, g.exlow, g.exhigh, [-1.0, 0.1], g.eyhigh, g.kind, g.label, g.draw, g.errors, g.style)
    @test_throws ArgumentError ROOTPlots.validate_graph(bad; label=false)
end

@testset "No request writer bypasses the Unicode translation (no ROOT)" begin
    # Every caller text reaches a request through root_text / scene_text (both
    # translate); clean_field alone is only used inside those two functions and in
    # the GraphSeries constructor, which checks but writes nothing.
    src = joinpath(pkgdir(ROOTPlots), "src")
    uses = String[]
    for file in readdir(src; join=true), (n, line) in enumerate(eachline(file))
        occursin("clean_field(", line) && !occursin("function clean_field", line) && push!(uses, "$(basename(file)):$(strip(line))")
    end
    @test sort(uses) == sort(["rootlatex.jl:return clean_field(root_latex(value))", "scene.jl:clean_field(value)",
                              "graphfit.jl:clean_field(s.label)"])
end

@testset "Unicode labels reach every ROOT text field translated (no ROOT)" begin
    mktempdir() do dir
        g = GraphSeries([1.0, 2.0], [1.0, 2.0]; yerrors=[0.1, 0.1], label="μ data")
        spec = GraphSpec(title="Δφ ≤ π", xlabel="σ_m [GeV/c²] → Δφ", ylabel="α β",
                         texts=[PlotText(0.2, 0.8, "χ²/ndf ≈ 1"; ndc=true)],
                         lines=[PlotLine(1.0, 1.0, 2.0, 2.0; label="λ line")],
                         xaxis=AxisStyle(change_labels=LabelChange(1; text="η=0")),
                         cms_label=(lumi="138 fb⁻¹", energy="13 TeV", additional_info="ℓ⁺ℓ⁻ channel"))
        path = joinpath(dir, "graph.tsv")
        ROOTPlots.write_graph_spec(path, [g], spec, ROOTPlots.graph_limits([g], spec))
        text = read(path, String)
        @test isascii(text) && check_ascii_request(path) === nothing
        for expected in ("#Delta#phi #leq #pi", "#sigma_m [GeV/c^{2}] #rightarrow #Delta#phi", "#alpha #beta",
                         "#mu data", "#chi^{2}/ndf #approx 1", "#lambda line", "#eta=0", "138 fb^{-1}",
                         "#font[12]{l}^{+}#font[12]{l}^{-} channel")
            @test occursin(expected, text)
        end
        # Histogram bin labels, legend entries, statistics-box header, comparison panels.
        h = Hist1D([0.0, 1.0, 2.0], [3.0, 4.0], [3.0, 4.0])
        hspec = GraphSpec(xaxis=AxisStyle(bin_labels=["α", "β"]), legend_spec=LegendSpec(entries=[LegendEntry(:series, 1, "γ hist")]))
        hpath = joinpath(dir, "hist.tsv")
        limits = ROOTPlots.hist_draw_limits(h, hspec, :e1, GraphStyle(), "δ")
        ROOTPlots.write_hist_draw_spec(hpath, h, hspec, limits, :e1, "δ", GraphStyle();
                                       request=ROOTPlots.stats_request(; stats=true, name="μ events"))
        htext = read(hpath, String)
        @test isascii(htext)
        @test all(s -> occursin(s, htext), ("#alpha", "#beta", "#gamma hist", "#delta", "#mu events"))
        cpath = ROOTPlots.write_root_spec(joinpath(dir, "cmp.tsv"), [comparison(h, [Component("τ MC", h, "#336699")])],
                                          [PanelSpec("m_{ℓℓ} [GeV]"; lumi="1 fb⁻¹")], (0.5, 1.5))
        ctext = read(cpath, String)
        @test isascii(ctext) && occursin("m_{#font[12]{l}#font[12]{l}} [GeV]", ctext) &&
              occursin("#tau MC", ctext) && occursin("1 fb^{-1}", ctext)
        # The request guard refuses non-ASCII text that bypassed translation.
        raw = joinpath(dir, "raw.tsv")
        write(raw, "ROOTPlotsGraph\t1\nLABELS\tσ\n")
        @test_throws ArgumentError check_ascii_request(raw)
        # Unsupported text is rejected before any ROOT process starts.
        @test_throws ArgumentError render_root([g], GraphSpec(xlabel="\uc9c8\ub7c9"), joinpath(dir, "korean"))
        @test_throws ArgumentError render_root(h, GraphSpec(), joinpath(dir, "korean-hist"); draw=:e1, label="\ud55c\uae00")
        @test !isfile(joinpath(dir, "korean.png")) && !isfile(joinpath(dir, "korean-hist.png"))
    end
end

if TEST_NATIVE_ROOT
    @testset "Unicode labels in native ROOT output (dump)" begin
        g = GraphSeries([1.0, 2.0, 3.0], [1.0, 3.0, 2.0]; yerrors=[0.2, 0.2, 0.2], label="μ data")
        spec = GraphSpec(title="Δφ ≤ π", xlabel="σ_m [GeV/c²] → Δφ", ylabel="α β",
                         texts=[PlotText(0.2, 0.8, "χ² ≈ 1"; ndc=true)],
                         cms_label=(lumi="138 fb⁻¹", energy="13 TeV"), canvas=CanvasStyle(opt_title=1))
        h = Hist1D([0.0, 1.0, 2.0], [3.0, 4.0], [3.0, 4.0])
        hspec = GraphSpec(xlabel="m_{ℓℓ} [GeV]", ylabel="Events",
                          xaxis=AxisStyle(bin_labels=["α", "β"], change_labels=LabelChange(1; text="η")))
        boxes = Ref{Any}(nothing)
        canvases = with_axis_dump() do dir
            render_root([g], spec, joinpath(dir, "graph"))
            boxes[] = with_stats_dump() do sdir
                render_root(h, hspec, joinpath(sdir, "hist"); draw=:e1, label="μ events", stats=true)
            end
        end
        @test length(canvases) == 2
        pad = first(framed(canvases[1]))
        @test pad.axes["x"].title == "#sigma_m [GeV/c^{2}] #rightarrow #Delta#phi"
        @test pad.axes["y"].title == "#alpha #beta"
        texts = [l.text for p in canvases[1].pads for l in p.latex]
        @test "#chi^{2} #approx 1" in texts
        @test any(t -> occursin("138 fb^{-1}", t), texts)
        @test any(p -> p.title[] !== nothing && p.title[].text == "#Delta#phi #leq #pi", canvases[1].pads)
        hpad = first(framed(canvases[2]))
        @test hpad.axes["x"].title == "m_{#font[12]{l}#font[12]{l}} [GeV]"
        @test hpad.axes["x"].bin_labels == Dict(1 => "#alpha", 2 => "#beta")
        @test only(hpad.axes["x"].changes).text == "#eta"
        @test only(boxes[]).lines[1] == "#mu events"
    end

    @testset "Mixed ROOT #-syntax and Unicode in native ROOT output (dump, widths)" begin
        # Unicode labels mixed with ROOT #-syntax reach ROOT exactly as expected, and
        # ROOT draws the commands: each translated label is narrower than the same
        # text with the command names spelled out and than the regressed "#^{}name".
        labels = ["#sigma_{m} [GeV/c²]", "Δ#phi", "x #in A, α", "#it{p}_{T} [μm]", "#splitline{α}{β}"]
        regressed = ["#^{}sigma_{m} [GeV/c^{2}]", "#Delta#^{}phi", "x #^{}in A, #alpha", "#^{}it{p}_{T} [#mum]",
                     "#^{}splitline{#alpha}{#beta}"]
        spelled = ["sigma_{m} [GeV/c^{2}]", "#Deltaphi", "x in A, #alpha", "it{p}_{T} [#mum]", "splitline{#alpha}{#beta}"]
        expected = Dict(MIXED_ROOT_LABELS)
        g = GraphSeries([1.0, 2.0, 3.0], [1.0, 3.0, 2.0]; yerrors=[0.2, 0.2, 0.2], label="#sigma data, σ")
        texts = PlotText[]
        for (k, set) in enumerate((labels, regressed, spelled)), (i, t) in enumerate(set)
            push!(texts, PlotText(0.15 + 0.27 * (k - 1), 0.85 - 0.12 * i, t; ndc=true, size=0.035))
        end
        push!(texts, PlotText(0.15, 0.15, "#alpha"; ndc=true, size=0.035))
        spec = GraphSpec(xlabel="#sigma_{m} [GeV/c²]", ylabel="Δ#phi", texts=texts, legend=:topright)
        canvases = with_axis_dump() do dir
            render_root([g], spec, joinpath(dir, "mixed"))
            @test filesize(joinpath(dir, "mixed.png")) > 0
        end
        pad = first(framed(only(canvases)))
        @test pad.axes["x"].title == "#sigma_{m} [GeV/c^{2}]" && pad.axes["y"].title == "#Delta#phi"
        @test any(e -> first(e) == "#sigma data, #sigma", only(pad.legends).entries)
        box = Dict(l.text => l for l in pad.latex)
        for i in eachindex(labels)
            fixed = box[expected[labels[i]]]
            @test fixed.width < box[spelled[i]].width
            @test fixed.width < box[regressed[i]].width
        end
        # #splitline draws two lines: the box is much taller than one Greek letter.
        @test box["#splitline{#alpha}{#beta}"].height > 1.5 * box["#alpha"].height
    end

    @testset "fit_root accepts Unicode graph labels (native)" begin
        x = collect(1.0:10)
        for label in ("σ data", "#sigma data, σ [GeV/c²]")
            fit = fit_root(GraphSeries(x, 2 .* x .+ 1; yerrors=fill(0.5, 10), label=label), :pol1; initial=[1.0, 1.0])
            @test fit.status == 0 && fit.valid
            @test fit.parameters ≈ [1.0, 2.0] atol=1e-6
        end
        # The fit ignores the presentation label: text TLatex cannot draw (Hangul,
        # emoji) is fitted unchanged, while drawing still rejects it clearly.
        reference = fit_root(GraphSeries(x, 2 .* x .+ 1; yerrors=fill(0.5, 10)), :pol1; initial=[1.0, 1.0])
        for label in ("\uce21\uc815 \ub370\uc774\ud130", "\U1F600 data", "#sigma \uac00")
            g = GraphSeries(x, 2 .* x .+ 1; yerrors=fill(0.5, 10), label=label)
            before = deepcopy(g)
            fit = fit_root(g, :pol1; initial=[1.0, 1.0])
            @test fit.status == 0 && fit.valid
            @test fit.parameters == reference.parameters && fit.covariance == reference.covariance
            @test g.label == label && g.x == before.x && g.y == before.y && g.eylow == before.eylow
            err = try render_root([g], GraphSpec(), joinpath(mktempdir(), "unsupported")); nothing catch e; e end
            @test err isa ArgumentError && occursin("ROOT TLatex cannot draw", err.msg)
        end
        # Request-breaking characters are still rejected when the series is built.
        @test_throws ArgumentError GraphSeries(x, x; label="a\tb")
        @test_throws ArgumentError GraphSeries(x, x; label="a\0b")
    end
end
