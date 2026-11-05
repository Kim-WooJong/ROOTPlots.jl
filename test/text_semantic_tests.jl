# Text correctness from ROOT's own state, independent of fonts and pixels.
# A one-character change of an axis title, a legend label, or the CMS
# luminosity text stays inside the pixel tolerance of the visual references, so
# the texts ROOT actually holds (ROOTPLOTS_AXIS_DUMP: axis titles, TLatex,
# legend entries; ROOTPLOTS_STATS_DUMP: statistics-box lines) are compared
# exactly. No geometry or font metrics are compared, so these checks run on
# every platform with native ROOT, also where the pixel references are skipped.

"""Normalized texts of a dumped figure: one line per pad with text, in drawing order."""
function text_snapshot(canvases)
    lines = String[]
    for canvas in canvases, pad in canvas.pads
        parts = String[]
        for a in ("x", "y", "z")
            haskey(pad.axes, a) && push!(parts, "$a=$(repr(pad.axes[a].title))")
        end
        pad.title[] === nothing || push!(parts, "title=$(repr(pad.title[].text))")
        append!(parts, ["latex=$(repr(t.text))" for t in pad.latex])
        append!(parts, ["legend=$(repr(e[1]))/$(e[2])" for l in pad.legends for e in l.entries])
        isempty(parts) || push!(lines, join(parts, " "))
    end
    return lines
end

const VISUAL_TEXTS = Dict(
    "th1_e1" => [raw"""x="x" y="Counts" """],
    "graph_band_tf1" => [raw"""x="x" y="y" legend="Data"/f legend="Model"/l """],
    "th2_colz" => [raw"""x="x" y="y" z="z" """],
    "thstack" => [raw"""x="x" y="Counts" legend="A"/f legend="B"/f """],
    "ratio_panel" => [raw"""x="" y="Events" latex="Synthetic" legend="A"/f legend="B"/f legend="Prediction stat. unc."/f legend="Data"/pe """,
                      raw"""x="x" y="Data / pred." """],
    "figure_2x2" => [raw"""x="x" y="y" """, raw"""x="x" y="y" z="z" """,
                     raw"""x="x" y="y" legend="A"/f legend="B"/f """,
                     raw"""x="x" y="y" legend="Data"/f legend="Model"/l """],
    "cms_style" => [raw"""x="x" y="Counts" latex="1 fb^{-1} (13 TeV)" latex="CMS" latex="Preliminary" """],
    "cms_outside" => [raw"""x="Integrated charge [C/cm^{2}]" y="Effective gain" latex="University of Seoul" latex="CMS" latex="Preliminary" legend="Measurement"/p legend="Mean = 6.14 #times 10^{4}"/l """],
)

@testset "texts of the visual reference figures (no pixels)" begin
    if TEST_NATIVE_ROOT
        cases = visual_cases()
        @test Set(first.(cases)) == Set(keys(VISUAL_TEXTS))   # a new visual case needs its texts here
        for (name, draw) in cases
            snapshot = with_axis_dump() do dir
                draw(joinpath(dir, name))
            end
            got = text_snapshot(snapshot)
            expected = rstrip.(VISUAL_TEXTS[name])
            got == expected || @warn "text snapshot differs" name got expected
            @test got == expected
        end
    else
        @test_skip false
    end
end

@testset "statistics-box texts and changed labels (no pixels)" begin
    if TEST_NATIVE_ROOT
        x = [0.5, 1.5, 1.5, 2.5, 2.5, 2.5, 3.5, -1.0, 9.0]
        h = histogram(x, 0:4).hist
        boxes = with_stats_dump() do dir
            render_root(h, GraphSpec(xlabel="x", ylabel="Events"), joinpath(dir, "stats"); draw=:hist,
                        stats=(:name, :entries, :underflow, :overflow), stats_name="Sample")
            render_root(merge(h, histogram([0.5], 0:4).hist), GraphSpec(), joinpath(dir, "merged"); draw=:hist,
                        stats=(:entries,))
        end
        # ROOT pads the box lines with spaces; whitespace runs are normalized to one space.
        normalized(box) = [join(split(l), " ") for l in box.lines]
        @test normalized(boxes[1]) == ["Sample", "Entries = 9", "Underflow = 1", "Overflow = 1"]
        @test normalized(boxes[2]) == ["Entries = 10"]
        # The changes that pass the pixel tolerance are caught here.
        base = GraphSpec(xlabel="x", ylabel="Counts", size=(500, 400), legend=:none, style=:CMS,
                         cms_label=(lumi="1 fb^{-1}", energy="13 TeV"))
        variants = (base, ROOTPlots.with_spec_fields(base; xlabel="y"), ROOTPlots.with_spec_fields(base; ylabel="Events"),
                    ROOTPlots.with_spec_fields(base; cms_label=(lumi="7 fb^{-1}", energy="13 TeV")))
        hv = Hist1D([0., 1., 2., 3., 4., 5.], [3., 7., 9., 6., 2.], [3., 7., 9., 6., 2.])
        snapshots = [text_snapshot(with_axis_dump() do dir
                         render_root(hv, spec, joinpath(dir, "v"); draw=:e1)
                     end) for spec in variants]
        @test allunique(snapshots)
        @test occursin("latex=\"7 fb^{-1} (13 TeV)\"", only(snapshots[4]))
        @test occursin("x=\"y\"", only(snapshots[2])) && occursin("y=\"Events\"", only(snapshots[3]))
    else
        @test_skip false
    end
end
