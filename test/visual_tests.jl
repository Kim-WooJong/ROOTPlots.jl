# Visual regression tests: render representative figures and compare them with
# reference PNGs in test/references/ using ROOT's own image reader.
#
# ROOTPLOTS_VISUAL selects the mode:
#   auto     (default) compare only when OS and ROOT version match the recorded
#            reference environment; otherwise the comparison is marked skipped.
#   required always compare; fail if references are missing or differ.
#   update   re-render and overwrite the references (review the diff before committing).
#   off      do not run.
# Fonts and anti-aliasing differ across platforms and ROOT releases, so pixel
# references are only meaningful in the environment that produced them.

const VISUAL_MODE = get(ENV, "ROOTPLOTS_VISUAL", "auto")
VISUAL_MODE in ("auto", "required", "update", "off") ||
    error("ROOTPLOTS_VISUAL must be auto, required, update, or off")
const REFERENCE_DIR = joinpath(@__DIR__, "references")
const REFERENCE_ENV = joinpath(REFERENCE_DIR, "environment.txt")
# Tolerance: anti-aliasing noise only. Real layout/color/marker changes exceed it.
const VISUAL_MAX_CHANGED_FRACTION = 0.002
const VISUAL_MAX_MEAN = 0.5

visual_environment() = "os=$(Sys.KERNEL)\nroot=$(strip(read(Cmd([ROOTPlots.find_root_config(), "--version"]), String)))\n"

"""Representative figures. Small canvases keep the suite fast; each covers a different painter path."""
function visual_cases()
    h = Hist1D([0., 1., 2., 3., 4., 5.], [3., 7., 9., 6., 2.], [3., 7., 9., 6., 2.])
    x = collect(0.:0.5:4.)
    data = GraphSeries(x, 1 .+ 1.5 .* x .+ 0.3 .* sin.(3 .* x); yerrors=fill(0.3, length(x)),
                       label="Data", errors=:band, style=GraphStyle(fill_color="#4477aa", fill_alpha=.3))
    model = ModelCurve(:pol1, [1., 1.5]; range=(0., 4.), label="Model")
    grid = Grid2D([0., 1., 2., 3.], [0., 1., 2., 3.], [1. 2. 3.; 2. 4. 6.; 3. 6. 9.])
    a = Hist1D([0., 1., 2., 3.], [4., 6., 3.], [4., 6., 3.])
    b = Hist1D([0., 1., 2., 3.], [2., 3., 5.], [2., 3., 5.])
    stack = HistStack([a, b]; labels=["A", "B"], fill_colors=["#cc6677", "#4477aa"])
    observed = Hist1D([0., 1., 2., 3.], [7., 8., 9.], [7., 8., 9.])
    cmp = comparison(observed, [Component("A", a, "#cc6677"), Component("B", b, "#4477aa")])
    small = GraphSpec(size=(500, 400), legend=:topleft)
    return [
        "th1_e1" => stem -> render_root(h, GraphSpec(xlabel="x", ylabel="Counts", size=(500, 400), legend=:none), stem; draw=:e1),
        "graph_band_tf1" => stem -> render_root(PlotPanel([data, model], small), stem),
        "th2_colz" => stem -> render_root(grid, Graph2DSpec(draw=:heatmap, size=(500, 400)), stem),
        "thstack" => stem -> render_root(stack, GraphSpec(xlabel="x", ylabel="Counts", size=(500, 400)), stem),
        "ratio_panel" => stem -> render_root([cmp], [PanelSpec("x"; annotation="Synthetic")], stem),
        "figure_2x2" => stem -> render_root(PlotFigure([PlotPanel(h, small), PlotPanel(grid, Graph2DSpec()),
                                                        PlotPanel(stack, small), PlotPanel([data, model], small)];
                                                       layout=(2, 2), size=(800, 600)), stem),
        "cms_style" => stem -> render_root(h, GraphSpec(xlabel="x", ylabel="Counts", size=(500, 400), legend=:none,
                                                        style=:CMS, cms_label=(lumi="1 fb^{-1}", energy="13 TeV")),
                                           stem; draw=:e1),
    ]
end

@testset "visual regression" begin
    if !TEST_NATIVE_ROOT || VISUAL_MODE == "off"
        @test_skip false
    else
        mktempdir() do dir
            # The comparator itself must detect real changes, independent of references.
            h = Hist1D([0., 1., 2.], [2., 3.], [2., 3.])
            spec = GraphSpec(size=(400, 300), legend=:none)
            base = render_root(h, spec, joinpath(dir, "base"); draw=:e1).png
            again = render_root(h, spec, joinpath(dir, "again"); draw=:e1).png
            moved = render_root(Hist1D([0., 1., 2.], [3., 2.], [3., 2.]), spec, joinpath(dir, "moved"); draw=:e1).png
            same = ROOTPlots.image_difference(base, again)
            @test same.same_size && same.changed_fraction <= VISUAL_MAX_CHANGED_FRACTION
            changed = ROOTPlots.image_difference(base, moved)
            @test changed.same_size && changed.changed_fraction > VISUAL_MAX_CHANGED_FRACTION

            recorded = isfile(REFERENCE_ENV) ? read(REFERENCE_ENV, String) : ""
            current = visual_environment()
            if VISUAL_MODE == "update"
                mkpath(REFERENCE_DIR)
                for (name, draw) in visual_cases()
                    cp(draw(joinpath(dir, name)).png, joinpath(REFERENCE_DIR, name * ".png"); force=true)
                end
                write(REFERENCE_ENV, current)
                @info "Visual references updated" REFERENCE_DIR
            elseif VISUAL_MODE == "auto" && recorded != current
                @info "Visual references were recorded elsewhere; comparison skipped" recorded current
                @test_skip false
            else
                for (name, draw) in visual_cases()
                    reference = joinpath(REFERENCE_DIR, name * ".png")
                    @test isfile(reference)
                    isfile(reference) || continue
                    candidate = draw(joinpath(dir, name)).png
                    d = ROOTPlots.image_difference(reference, candidate)
                    ok = d.same_size && d.changed_fraction <= VISUAL_MAX_CHANGED_FRACTION &&
                         d.mean <= VISUAL_MAX_MEAN
                    ok || @warn "Visual difference" name d candidate
                    @test ok
                end
            end
        end
    end
end
