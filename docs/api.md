# ROOTPlots API and computational semantics

[Quick start](../README.md) · [Installation](installation.md) · [Migration](migration.md) · [Testing](testing.md)

This document describes the supported scope and the scientific contracts of the
computations. ROOTPlots is not a complete ROOT API binding.

ROOTPlots.jl is a Julia API for using ROOT's own drawing and fitting features,
not a reimplementation. Every option maps onto a ROOT method/option/style
setting, named so a ROOT user recognises it (snake_case of the ROOT method;
docstring cites it). An unset option keeps ROOT's default behaviour. Drawing,
statistics, fitting and formatting are done by ROOT; Julia only validates input
and passes data/options. See [design](design.md).

## Supported plots

| Task | Julia API | Native ROOT object / behavior |
|---|---|---|
| 1D histograms: errors, bars, steps, overlays | `Hist1D`, `histogram`, `render_root` | TH1D, ROOT draw options |
| Weighted 2D/3D histograms with errors | `Hist2D`, `Hist3D` | TH2D/TH3D contents and Sumw2 errors |
| Heatmaps, contours, LEGO, surfaces, volumes | `Grid2D/3D`, `Graph2DSpec/Graph3DSpec` | TH2D/TH3D painters |
| Stacked / unstacked histograms | `HistStack` | THStack; data/model comparison panels |
| 1D/2D/3D profiles | `Profile1D/2D/3D`, `profile_stats` | TProfile/TProfile2D/TProfile3D filled natively |
| Polygon bins | `PolygonHist2D` | TH2Poly (simple convex or non-convex polygons, no holes) |
| XY points, lines, curves, fills, bars | `GraphSeries`, `GraphStyle` | TGraph and native draw options |
| Symmetric/asymmetric errors, bands, boxes, arrows | `GraphSeries(...; xerrors, yerrors, errors)` | TGraphErrors/TGraphAsymmErrors |
| Bent (shifted) error bars | `BentErrorSeries` | TGraphBentErrors |
| Several y-error sets per point (e.g. stat + syst) | `MultiErrorSeries`, `YErrorSet`, `total_yerrors` | TGraphMultiErrors |
| XYZ points/lines, errors, explicit interpolation | `Graph3D`, `Graph3DErrors` | TGraph2D/TGraph2DErrors |
| Sampled Julia functions / native TF1 | `function_curve` / `ModelCurve` | Sampled TGraph / TF1 of a predefined model |
| Native TF2 / TF3, sampled Julia functions of 2–3 variables | `ModelSurface` / `ModelVolume` / `function_grid` | TF2, TF3 of predefined models; sampled TH2D/TH3D |
| Polar, QQ, spline | `PolarGraph`, `QQGraph`, `SplineCurve` | TGraphPolar, TGraphQQ, TSpline3/5 |
| Pie charts, binomial efficiency | `PieChart`, `EfficiencyGraph` | TPie, TEfficiency with native intervals |
| Lines, LaTeX text, legends, axes, margins | `PlotLine`, `PlotText`, `LegendSpec`, `AxisStyle`, `CanvasStyle` | TLine/TLatex/TLegend/ROOT axes |
| Arrows, boxes, ellipses, polygons, secondary axes | `PlotArrow`, `PlotBox`, `PlotEllipse`, `PlotPolygon`, `SecondaryAxis` | TArrow/TBox/TEllipse/TPolyLine/TGaxis |
| Mixed and nested multi-panel figures, vector PDF | `PlotPanel`, `PlotFigure` | One TCanvas with (nested) TPads (not stitched PNGs) |
| Parallel work | `rootmap`; `histogram(...; threads=...)` | Julia tasks/CPU threads + isolated ROOT processes |
| Statistics and fit boxes | `stats`, `fitstats`, `attached_fit`, `stats_position` keywords | TPaveStats from `SetOptStat`/`SetOptFit`, TF1 with fit results |

Out of scope or not implemented: ROOT file I/O, RDataFrame, GUI/OpenGL,
RooFit/RooStats, arbitrary C++ or TF1/TF2/TF3 expression strings (analytic
functions come from fixed named formulas), OpenGL-only TF3/TH3 views, polygons
with holes, and geometry/event displays. Not every axis detail of polar/XYZ
plots and not every ROOT class method is exposed. A Julia function can be
sampled and overlaid, but that is not described as the corresponding special
ROOT object.

The supported target is ROOT 6.40 or newer; no compatibility layer for older
ROOT 6 releases is provided. See [testing](testing.md) for the environments that
have actually been exercised.

### Histograms, profiles, stacks

```julia
using ROOTPlots

one = histogram([0.2, 0.4, 1.2], [0., 1., 2.]; weights=[1., 2., 1.]).hist
two = histogram(([0.2, 0.4, 1.2], [0.3, 1.3, 1.4]), ([0., 1., 2.], [0., 1., 2.]);
                weights=[1., 2., 1.], threads=:auto)
render_root(two.hist, Graph2DSpec(draw=:heatmap), "output/h2"; format=:both)
px = projection(two.hist, :x)  # sums visible bins; sumw2 is summed separately

stack = HistStack([one, one]; labels=["A", "B"], fill_colors=["#cc6677", "#4477aa"])
render_root(stack, GraphSpec(xlabel="x", ylabel="Counts"), "output/stack")

profile = Profile1D([0., 1., 2.], [.2, .4, 1.2, 1.4], [1., 3., 2., 5.]; error=:mean)
render_root(profile, GraphSpec(xlabel="x", ylabel="Mean y"), "output/profile")
statistics = profile_stats(profile)  # visible-bin mean/error/entries/effective_entries
```

Inputs may be 1D vectors, DataFrame columns, 2D/3D coordinate tuples, or tuples
of DataFrame column names. For N-dimensional results, `bins_with_flow=(sumw, sumw2)`
holds the **complete** arrays including underflow/overflow on every axis, while
`hist` holds only visible bins. As in ROOT, the upper edge belongs to the
overflow. Threaded filling accumulates privately per chunk and merges in a fixed
order, preserving signed weights. `projection(h, (:z, :x))` orders the
resulting 2D axes as requested. No automatic unit-area or density normalization
is applied.

Profiles own their raw coordinates, measurements, and optional `weights` and
are computed with ROOT's `Fill`. `error=:mean|:spread|:integer|:gaussian` is
supported with finite non-negative weights; Gaussian mode requires positive
weights. A true stack requires identical binning and non-negative contents;
signed bins are only accepted with `draw="NOSTACK HIST"`. `PolygonHist2D`
displays pre-aggregated values for simple polygons (convex or not; see
[polygons](#non-convex-polygons)) and does not fill events.

### Mixed panels and annotations

```julia
g = GraphSeries([0., 1., 2.], [1., 2., 3.]; yerrors=[.1, .2, .2], label="Data",
                errors=:band, style=GraphStyle(fill_color="#4477aa", fill_alpha=.3))
model = ModelCurve(:pol1, [1., 1.]; range=(0., 2.), label="TF1")
figure = PlotFigure([
    PlotPanel([g, model], GraphSpec(xlabel="x", ylabel="Response");
              primitives=[PlotArrow(.3, 2.5, 1., 2.)]),
    PlotPanel(two.hist, Graph2DSpec(title="Weighted TH2");
              zlimits=(0., 4.), zaxis=AxisStyle(title_size=.06), contour_levels=[1., 2., 3.]),
    PlotPanel(stack, GraphSpec(title="THStack")),
    PlotPanel(profile, GraphSpec(title="TProfile"))
]; layout=(2, 2), size=(1200, 900))
render_root(figure, "output/figure"; format=:both)
```

Panels are row-major; `PlotFigure([p11 p12; p21 p22])` is also accepted. XY
panels may mix `Hist1D`, `GraphSeries`, `BentErrorSeries`, `MultiErrorSeries`, and
`ModelCurve` layers; ND panels
draw one ND object. The `PlotPanel` keywords `xaxis/yaxis/zaxis`,
`xlimits/ylimits/zlimits`, `canvas`, `lines/texts/legend_spec`,
`style/styles`, and `label/labels` apply to basic XY/ND panels. Advanced,
ModelCurve, TF2/TF3, XYZ, and special objects use their own spec and style; common
`primitives` may be added to any panel. Because a basic XY panel owns its axes,
use standalone `render_root` for special TH1 options such as HBAR/PIE/FUNC/AXIS.

Palettes and error-cap/PostScript line scale are ROOT global settings, so
conflicting values within one figure are rejected; ND panels must share a
palette (this includes panels in nested sub-grids). Automatic shared axes are
not provided. For grid panels, `zlimits` is the color range for TH2 and the
spatial z range for TH3.

#### Nested layouts

A cell of a `PlotFigure` may itself be a `PlotFigure`; it is drawn as a
sub-grid inside that cell (ROOT `TPad::Divide` of the cell's pad, i.e. genuinely
nested TPads). `widths`/`heights` give relative column/row weights.

```julia
right = PlotFigure([PlotPanel(two.hist, Graph2DSpec(title="TH2")),
                    PlotPanel(profile, GraphSpec(title="TProfile"))]; layout=(2, 1), heights=[1, 2])
figure = PlotFigure([PlotPanel([g, model], GraphSpec()), right];
                    layout=(1, 2), widths=[3, 2], size=(1300, 700))
render_root(figure, "output/nested")
```

Figures nest at most 4 levels deep and hold at most 64 panels in total, with
1:16 rows/columns per grid and weight ratios up to 100. Only the outermost
figure's `size` is used. Without `widths`/`heights` and nested cells, the
native request is identical to the previous flat layout. ROOT text and marker
sizes are fractions of the pad, so labels in small sub-pads become small;
enlarge them with `AxisStyle`/`PlotText` sizes if needed.
Arrow/Box/Ellipse/Polygon/SecondaryAxis use data coordinates; a secondary axis is a
linear numeric mapping and does not transform data.

### Functions and special graphs

```julia
sampled = function_curve(x -> sin(x), (0., 2pi); n=400)        # Julia callable → TGraph
analytic = ModelCurve(:gaus, [1., 0., 1.]; range=(-3., 3.), npx=400)  # native TF1
spline = SplineCurve([0., 1., 2., 3.], [0., 1., .5, 2.]; order=3, npx=400)
efficiency = EfficiencyGraph([0., 1., 2.], [5., 9.], [10., 10.])
render_root(efficiency, GraphSpec(xlabel="x", ylabel="Efficiency"), "output/efficiency")
points = efficiency_points(efficiency)  # GraphSeries of ROOT Clopper–Pearson intervals
```

`errors=:bars|:none|:band|:smoothband|:boxes|:caps|:nocaps|:arrows|:filledarrows`
and `draw=:points|:line|:linepoints|:curve|:curvepoints|:fill|:bar|:bar_offset`
are supported. Raw graph option strings may not contain `A` or `SAME` (the
frame owns the axes; `SAME` is appended internally), and a raw string takes
precedence over `errors`. Filled or smoothed bands are a presentation choice,
not a new confidence-interval computation.

`PolarGraph` angles are in radians; Cartesian log/range settings are rejected.
`QQGraph(x, y)` puts the first equal-length sample on the horizontal axis and
the second on the vertical axis (the ROOT C++ constructor order is reversed
internally); `qq_points` returns ROOT's quantile points. `SplineCurve`
requires sorted distinct x and supports cubic or quintic splines; `spline_points`
evaluates the native spline for display. `EfficiencyGraph` uses unweighted
integer `passed <= total` with Clopper–Pearson intervals (configurable
`confidence`); bins with zero denominator are omitted. `PieChart(values; labels,
colors)` is a native TPie.

### Bent and multiple error bars

```julia
x = collect(1.:6.)
a = BentErrorSeries(x, ya; yerrors=(low=fill(.45, 6), high=fill(.35, 6)), ybend=fill(-.18, 6), label="A")
b = BentErrorSeries(x, yb; yerrors=fill(.4, 6), ybend=fill(.18, 6), label="B")
render_root([a, b], GraphSpec(xlabel="x", ylabel="y"), "output/bent")

stat = YErrorSet(fill(.25, 6); label="Statistical")
syst = YErrorSet((low=sl, high=sh); label="Systematic", display=:boxes, xscale=.6,
                 style=GraphStyle(color="#4477aa", fill_alpha=.45))
data = MultiErrorSeries(x, ya, [stat, syst]; xerrors=fill(.3, 6), label="Data", sum=:quadrature)
render_root(data, GraphSpec(xlabel="x", ylabel="y"), "output/stat-syst")
combined = total_yerrors(data)   # (low=…, high=…) with the selected sum rule
```

`BentErrorSeries` is a native `TGraphBentErrors`: the error magnitudes are as for
`GraphSeries`, and `ybend`/`xbend` move the error-bar **end points** (in axis
units) without moving the point: the lower y bar ends at `(x + ybend.low, y - eylow)`,
the left x bar at `(x - exlow, y + xbend.low)`. Use it to separate overlapping
bars; the bend is display geometry only. `errors` accepts the `GraphSeries` modes.

`MultiErrorSeries` is a native `TGraphMultiErrors` with 1:16 y-error sets
(`YErrorSet`), each drawn separately from the point with its own line/fill style
and `display` (`:bars`, `:nocaps`, `:caps`, `:arrows`, `:filledarrows`, `:boxes`,
`:boxes_outline`, `:band`, `:smoothband`, `:none`). Box widths are the x
errors times `xscale`. Labelled sets get their own legend entries after the
regular layers (explicit `LegendEntry(:series, i, …)` indices are unchanged).
`sum=:first|:quadrature|:linear` is ROOT's `SetSumErrorsMode`
(`kOnlyFirst/kSquareSum/kAbsSum`); it does not change how the sets are drawn and
is reproduced by `total_yerrors`. Whether a quadrature sum is justified
(independent sources) is the caller's decision. Neither type is a fit input;
fit a `GraphSeries` built with the chosen combined errors instead.

Both types work as standalone figures, in a vector of mixed XY layers
(`render_root([graph, bent, multi, model], spec, stem)` draws them on one pad),
and in `PlotPanel`. Frame ranges include all drawn end points (bent offsets,
the widest visible error set); on log axes these end points must be positive.

### Analytic TF2/TF3 and sampled Julia functions

Like `ModelCurve` (TF1), the 2-D/3-D analytic objects are genuine ROOT functions
built from a **fixed list of named formulas** defined on both the Julia and C++
side; no user expression string reaches TFormula/Cling. Arbitrary Julia code
is supported by sampling instead.

```julia
surface = ModelSurface(:bigaus, [1., 0., 1., 0., .7, .5]; xrange=(-3, 3), yrange=(-3, 3), npx=100, npy=100)
render_root(surface, Graph3DSpec(draw="SURF1"), "output/tf2-surface")
render_root(surface, Graph2DSpec(draw="CONT1Z", contours=12), "output/tf2-contours")
render_root(surface, Graph2DSpec(draw=:heatmap), "output/tf2-colz")

volume = ModelVolume(:gaus3d, [1., 0., 1., 0., .7, 0., 1.3];
                     xrange=(-3, 3), yrange=(-3, 3), zrange=(-3, 3), level=.3)
render_root(volume, Graph3DSpec(draw=:iso), "output/tf3-iso")   # surface f = 0.3
render_root(volume, Graph3DSpec(draw=:box), "output/tf3-box")   # sampled TH3D, f > 0.3

grid = function_grid((x, y) -> sin(x) * cos(y), (0, pi), (0, pi); n=(100, 100))  # Grid2D
```

| Model | Parameters (1-based) |
|---|---|
| `:gaus2d` | `[A, mx, sx, my, sy]`: `A exp(-((x-mx)/sx)²/2 - ((y-my)/sy)²/2)` |
| `:bigaus` | `[A, mx, sx, my, sy, rho]`, correlated Gaussian, `abs(rho) < 1` |
| `:pol1_2d` / `:pol2_2d` | `[p0, px, py]` / `[p0, px, py, pxx, pxy, pyy]` |
| `:gaus3d` | `[A, mx, sx, my, sy, mz, sz]` |
| `:pol1_3d` / `:pol2_3d` | `[p0, px, py, pz]` / `[p0, px, py, pz, pxx, pyy, pzz, pxy, pxz, pyz]` |

`ModelSurface` (TF2) accepts `Graph2DSpec` (`:heatmap` → `COLZ`, `:contour` →
`CONT3`, or `COL[0-2][Z]`, `CONT[0-4][Z]`) and `Graph3DSpec` (`:surface` →
`SURF2Z`, or `SURF[1-7][Z]`, `LEGO[1-4][Z]`); log axes need positive ranges
(and positive values for `logz`). ROOT paints a TF2 by evaluating it on an
`npx × npy` grid, so contours and surfaces are grid renderings of the exact
function. `ModelVolume` (TF3) offers only what ROOT draws in batch mode
without OpenGL: `:iso`, ROOT's implicit-surface painter for `f(x, y, z) = level`,
and `:box`, the TF3 evaluated natively at the centres of an `npx × npy × npz`
TH3D drawn as `BOX2Z` with bins `f ≤ level` hidden. GL options, transparent
volume rendering, slices/projections of TF3, log axes for TF3, and fitting with
TF2/TF3 are not provided.

`function_grid(f, xrange, yrange; n)` / `function_grid(f, xrange, yrange, zrange; n)`
evaluate a Julia callable at the cell centres and return a `Grid2D`/`Grid3D`
(TH2D/TH3D), the grid analogue of `function_curve`. All `Grid2D`/`Grid3D`
options apply; this is a display-only sampling, not a TF2/TF3.

### Non-convex polygons

`PolygonHist2D` bins and the `PlotPolygon` primitive accept any **simple**
polygon, convex or not (L shapes, stars, U shapes):

```julia
L = ([0., 2, 2, 1, 1, 0], [0., 0, 1, 1, 3, 3])
cells = PolygonHist2D([L, star], [1., 3.])
render_root(cells, Graph2DSpec(), "output/polygons")
panel = PlotPanel([graph], GraphSpec(); primitives=[PlotPolygon(L...; fill_color="#cc6677", fill_alpha=.5)])
```

Each polygon needs 3:256 finite vertices (a repeated closing vertex is
removed), no repeated neighbours, a non-zero area, and no self-intersection:
non-adjacent edges may not touch or cross and an edge may not fold back onto
its neighbour. Violations raise `ArgumentError` naming the problem, and the
native bridge repeats the check. Holes (polygons with inner rings) and
multi-part bins are not supported. `PlotPolygon` is drawn as a filled
`TPolyLine` plus outline in data coordinates; `fill_style=0` draws the outline only.

Runnable examples: `scripts/graphics_gallery.jl` (9 panels),
`scripts/special_graph_demo.jl` (6 panels), `scripts/advanced_histogram_demo.jl`,
`scripts/graph_extension_demo.jl` (bent/multi errors, TF2/TF3, polygons, nested
pads), `scripts/edge_profile_demo.jl`. All inputs are synthetic.

## Parallel ROOT jobs

```sh
julia --threads=4 --project=. scripts/parallel_demo.jl
```

```julia
results = rootmap(eachindex(datasets); workers=4) do i
    data = datasets[i]
    fit = fit_root(data, :pol1; initial=[0., 1.])
    fit.status == 0 && fit.valid || error("fit $i failed")
    files = render_root([fit_curve(fit), data], GraphSpec(xlabel="x", ylabel="y"),
                        "output/measurement-$i"; format=:both)
    (; fit, files)
end
```

- `rootmap(f, jobs; workers=Threads.nthreads(:default))` runs at most `workers`
  Julia tasks on the default thread pool and returns results in **input order**.
  `workers=1` runs sequentially.
- **Event filling:** `histogram(...; threads=:auto|:on|:off)` is Julia CPU
  multithreading with private accumulators and a fixed merge order.
- **Fits and figures:** every ROOT call runs in its own C++ ROOT process, so ROOT
  global style/canvas state is isolated per process. A single canvas or fit is
  not parallelized by ROOT-internal threads.
- With one Julia thread, `workers > 1` still overlaps waiting on external
  processes. Start Julia with `--threads=N` to parallelize CPU-side preparation.
  `workers` is the number of concurrent jobs, not cores or ROOT threads.
- Give **every job a distinct output path**, and do not mutate shared data, spec
  arrays, DataFrames, or the process environment while jobs run. No deep copies or
  duplicate-path checks are performed; this also applies to concurrent calls
  outside a batch.
- All inputs are attempted and workers are joined before returning. Failures are
  thrown as `RootMapError` whose `failures` hold `index => CapturedException` in
  input order. Outputs completed earlier remain; a batch is not transactional.
  Callers must still check each fit's `status/valid`.
- Small jobs may not speed up because of process start-up cost. Limit `workers`
  according to per-job CPU/RAM use and nested histogram threading. No fixed
  speed-up is promised.

## Binning and data/model comparisons

```julia
using DataFrames, ROOTPlots

events = DataFrame(mass=[25.0, 45.0, 100.0], weight=[1.0, -0.2, 2.0])
binned = histogram(events, :mass, [20.0, 40.0, 60.0, 120.0]; weight=:weight, threads=:auto)
# binned.hist: visible bins only; binned.underflow / overflow: (Σw, Σw²) each
```

`histogram` reads DataFrame/SubDataFrame columns directly. Data and each model
component can be binned identically and combined in a comparison panel:

```julia
edges = [20.0, 40.0, 60.0, 120.0]
observed = histogram(data_events, :mass, edges).hist
model = histogram(model_events, :mass, edges; weight=:weight, threads=:auto).hist
cmp = comparison(observed, [Component("Model", model, "#C5001A")])
spec = PanelSpec("Mass [GeV]"; logy=true, ymin=0.1)
render_root([cmp], [spec], "output/mass")
```

The default output is `format=:png`. Use `format=:pdf` for PDF only and
`format=:both` for both; this applies to every `render_root` figure type. The
returned `png`/`pdf` keys are always present and the unselected one is `nothing`.
Existing files of an unselected format are **not deleted**. The native executable
writes both formats to a temporary directory and only the selected files are
moved to the destination.

When Julia is started with `--threads=auto` or `JULIA_NUM_THREADS=4`,
`threads=:auto` uses parallel filling for large enough inputs. `:on` requests
parallel execution where possible and `:off` processes the same chunks serially.
The keyword does not create threads. Because chunk accumulators merge in a
fixed order, all three modes are bit-identical in one environment; results may
differ by floating-point rounding from naive row-order accumulation. Do not
modify the DataFrame or referenced columns while filling.

Visible bins are `[low, high)` as in ROOT, and **the top edge is overflow**.
Whether to fold flow bins into the first/last visible bin is an explicit caller
decision. Components with different binning are never rebinned implicitly.

To compare histograms with different fine binning, merge them **explicitly**.
`rebin(h, new_edges)` requires the new edges to be an **exact subset** of the
existing edges with identical outer edges. `Σw` and `Σw²` are summed per new
bin, preserving negative weights. Bins are never split or cropped, and visible
content is never moved into flow bins. The returned `Hist1D` shares no arrays
with its input.

```julia
common = [0.0, 3.0, 6.0]
data = Hist1D([0.0, 1.0, 3.0, 6.0], [2.0, 1.0, 4.0], [2.0, 1.0, 4.0])
model = Hist1D([0.0, 2.0, 3.0, 6.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
cmp = comparison(rebin(data, common), [Component("Model", rebin(model, common), "#336699")])
```

`Hist1D` does not record whether contents are yields or densities. `rebin` merges
yields; summing densities would be wrong. Unit-area normalization by the data's
own total introduces bin-to-bin covariance and is intentionally not provided.

The stacked prediction is the per-bin sum of component `Σw`. Its default
statistical uncertainty is the square root of the summed component `Σw²`
(independent components). Passing `prediction_errors=BinErrors(...)` **replaces**
the displayed total uncertainty. Systematic uncertainties, correlations, and
post-fit covariances must be computed by the caller.

For `P > 0` the ratio is `D/P`, the data error is `data error / P`, and the
prediction band is `1 ± prediction error / P`. For `P ≤ 0` the ratio is undefined
and neither point nor band is drawn for that bin. `D = 0, P > 0` is valid;
asymmetric Poisson errors can be supplied via `BinErrors`. Negative **event
weights** are preserved during filling, but negative **component bin sums**
cannot be drawn as a conventional stack and are rejected.

## Graphs and fitting

```julia
x = collect(0.0:4.0)
y = [1.1, 2.9, 5.2, 6.8, 9.1]
observed = GraphSeries(x, y; yerrors=fill(0.2, 5), label="Measurement")
fit = fit_root(observed, :pol1; initial=[0.0, 1.0])
@show fit.parameters fit.errors fit.covariance fit.status fit.chi2 fit.ndf

curve = fit_curve(fit; n=300, label="ROOT fit")
render_root([curve, observed], GraphSpec(xlabel="x", ylabel="Response", legend=:topleft),
            "output/graph-fit")
```

`GraphSeries` maps to `TGraph`, `TGraphErrors`, or `TGraphAsymmErrors`.
`xerrors`/`yerrors` accept a symmetric vector, `(low=…, high=…)`, or `BinErrors`;
errors are non-negative **magnitudes**. Connected lines follow input order and
are never sorted or interpolated automatically. `xlimits`/`ylimits` in `GraphSpec`
are physical coordinates, and `logx`/`logy` are independent. On log axes, the
centers and the displayed error endpoints must be positive. A standalone call
writes a single pad; use `PlotPanel/PlotFigure` for multiple pads.

`fit_root` runs ROOT's `TH1::Fit`/`TGraph::Fit` and returns a `FitResult` with
parameters, standard errors, the full covariance, status, validity, `chi2`,
`ndf`, the number of points used, and the minimizer. Models are `:pol0`,
`:pol1`, `:pol2`, `:gaus` (`A exp(-0.5((x-μ)/σ)^2)`), `:exp` (`A exp(bx)`), and
the two edge models below. `initial` is required and follows ROOT's parameter
order. `bounds=[nothing, (low, high), ...]`, `fixed=[1, ...]`, and
`range=(low, high)` are optional. `fit_curve` samples the returned parameters
for display; it is not a confidence band.

| Model | Parameter order (1-based) | Function |
|---|---|---|
| `:erfc_edge` | `[B, A, R, sigma]` | `B + A/2 * erfc((x-R)/(sqrt(2)*sigma))` |
| `:double_erf` | `[B, H, xL, xR, sigmaL, sigmaR]` | `B + H/2 * (erf((x-xL)/(sqrt(2)*sigmaL)) - erf((x-xR)/(sqrt(2)*sigmaR)))` |

Widths must be positive. A free width requires explicit finite positive bounds
(for example `(0.03, 1.5)`); a width fixed to a positive value needs none. These
models are defined explicitly on both the Julia and C++ side instead of executing
arbitrary TF1 strings. `H` is a coefficient; when the edges overlap or are
asymmetric it may differ from the background-subtracted peak height.
`scripts/edge_profile_demo.jl` shows both models on synthetic data.

### Axes, margins, reference lines, LaTeX annotations, legends

```julia
observed = GraphSeries([0., 1., 2.], [8., 5., 2.]; yerrors=[1., 0.7, 0.4], label="Profile")
spec = GraphSpec(xlabel="Radius r (mm)", ylabel="Density (counts / mm^{2})",
    xlimits=(0., 2.5), ylimits=(0., 12.), size=(900, 720),
    xaxis=AxisStyle(title_font=42, label_font=42, title_size=0.06,
                    label_size=0.05, title_offset=1., divisions=505, max_digits=3),
    canvas=CanvasStyle(left_margin=0.14, right_margin=0.04,
                       bottom_margin=0.14, top_margin=0.05,
                       ticks_x=true, ticks_y=true, end_error_size=0.),
    lines=[PlotLine(1.5, 0., 1.5, 9.; color="#666666", style=2)],
    texts=[PlotText(0.58, 0.84, "#sigma = 180 #pm 5 #mum"; size=0.038)],
    legend_spec=LegendSpec(bounds=(0.18, 0.74, 0.52, 0.90), text_size=0.037,
        entries=[LegendEntry(:series, 1, "Profile"; option="lep"),
                 LegendEntry(:line, 1, "Nominal radius"; option="l")]))
render_root([observed], spec, "output/profile"; format=:both)
```

- `AxisStyle`: per-axis `title_font`, `title_size`, `title_offset`, `label_font`,
  `label_size`, `label_offset`, `tick_length`, `divisions`, `max_digits`,
  `center_title`, plus the number/label options of
  [Axes, titles and number format](#axes-titles-and-number-format). Omitted
  fields keep ROOT defaults; precision-3 fonts (e.g. 43) use pixel sizes.
- `CanvasStyle`: four margins, ticks on both sides, `frame_line_width`,
  `end_error_size`, `line_scale_ps`, grids, and the opt-in ROOT title box.
  Margins are canvas fractions.
- `PlotLine` uses **data coordinates** by default; `ndc=true` uses 0–1 canvas
  coordinates. A line `label` is collected into the automatic legend.
- `PlotText` uses **NDC coordinates** by default; `ndc=false` uses data
  coordinates. `font`, `size`, `color`, `align`, and `angle` are configurable.
  Math uses ROOT TLatex syntax (`#sigma`, `#pm`, `#mum`, `mm^{2}`), not a full
  LaTeX engine.
- An automatic legend (`legend=:topright` etc.) uses a fixed text size of
  0.035 of the pad (ROOT's default axis label size) and a box sized from its
  entries, so its text matches the axis labels in single canvases and figure
  panels alike.
- `legend_spec` replaces the `legend=:topright`-style position. `background=nothing`
  is transparent; a color string is opaque. `entries=nothing` collects labels
  automatically and `entries=[]` shows none. `LegendEntry` sets order, target, and
  the `l/p/pe/ep/lp/lep/f` marker.
- Annotations never extend the data range. Data-coordinate lines/text must be
  positive on log axes. Clipping of out-of-range annotations follows ROOT pads.
- Lines, texts, and legends apply to **single-pad** plots. Axis and canvas
  styles also apply to 2D/3D specs (`Graph2DSpec`/`Graph3DSpec` `xaxis`, `yaxis`,
  `zaxis`, `canvas`), comparison panels (`PanelSpec` `xaxis`, `yaxis`,
  `ratio_yaxis`, grids/ticks), and the two-pad fit/pull figure (axis styles, grids,
  ticks, title box); settings a layout owns (e.g. its margins) are rejected
  instead of silently ignored.

A native `TH1D` can share a canvas with graphs or fit curves:

```julia
h = Hist1D([0., 1., 2., 3.], [2., 3., 1.], [2., 3., 1.])
curve = GraphSeries([0., 1., 2., 3.], [1., 2.5, 2.5, 1.]; draw=:line, label="Model")
render_root(h, [curve], GraphSpec(xlabel="x", ylabel="Counts / bin"),
            "output/hist-and-model"; draw="E1", format=:both)
```

Here the legend `:series` index is 1 for the histogram and 2, 3, … for the
graphs. Automatic ranges include all layers. Fit histogram data as the original
`Hist1D` to keep TH1 χ²/likelihood semantics.

Histograms can be fitted directly:

```julia
h = Hist1D([0., 1, 2, 3, 4], [2., 0, 4, 2], [2., 0, 4, 2])
poisson_fit = fit_root(h, :pol0; initial=[1.0], objective=:poisson)
```

With the default `objective=:chi2`, ROOT may exclude zero-error bins. `:poisson`
accepts only unweighted integer counts (`sumw2 == sumw`); `:weighted_poisson` uses
weights `sumw/sumw2`. Both likelihood modes reject negative bins, and the weighted
mode requires positive `sumw2` for positive contents. Graph fits support `:chi2`
only, and supplied y errors must be strictly positive. ROOT's effective-variance
treatment of x errors is the default `xerrors=:effective`; use `:ignore` to
disable it. `chi2_kind=:poisson_deviance` is ROOT's likelihood deviance for
unweighted Poisson fits; weighted Poisson reports `:weighted_poisson_objective`.
**Do not judge a fit by `status == 0` alone**; check `valid`, the covariance
status, residuals, and model adequacy.

## Histogram statistics and fit diagnostics

```julia
h = histogram(events, :mass, edges; weight=:weight)
s = histstats(h.hist; underflow=h.underflow, overflow=h.overflow)
@show s.integral s.effective_entries s.mean s.mean_error s.std

fit = fit_root(h.hist, :pol0; initial=[1.0])
r = fit_residuals(h.hist, fit)
@show r.centers r.residual r.pull r.valid
@show fit_correlation(fit) chi2_pvalue(fit)
```

`integral(hist; bins=…, width=false)` returns `(value, error)` for the selected
**visible** bins; `width=true` multiplies contents and errors by bin width (use it
for **densities**). `effective_entries(hist)` is `(Σw)²/Σw²`, not the event count or
ROOT's `GetEntries()`. `binned_mean`, `binned_std`, and `binned_mean_error` use bin
centers; the last is ROOT's approximation `std/√effective_entries`. Exact event
moments cannot be recovered from bin contents alone (`histogram` records them in
`hist.moments` for ROOT's statistics box). `histstats` bundles these values and
reports underflow/overflow **separately without using them**. A non-positive
total weight, or a negative variance from signed weights, raises an error.

`fit_residuals` accepts χ² histogram fits only. It computes `observed − model` at
bin centers and the pull divided by `√sumw2`; bins outside the fit range or with
zero error have `valid=false` and `NaN` residual/pull. `fit_correlation` is the
correlation matrix; rows/columns of fixed or zero-error parameters are `NaN`.
`chi2_pvalue` returns the asymptotic upper-tail probability for valid χ² or
**unweighted** Poisson fits, and refuses weighted Poisson, zero ndf, or fits
with ROOT-normalized errors. Interpret it only under the usual assumptions:
sufficient counts, a fit range chosen beforehand, no parameter at a bound, and an
adequate model.

`histogram_points(h)` converts bins to **center points** with `√sumw2` vertical
errors for overlays. It is not a step histogram and bin widths are not used as x
errors. **Fit the original `Hist1D`** to keep TH1 semantics. `fit_pull_points(h, fit)`
plots the pulls of the bins used by a valid χ² histogram fit.

`histogram_steps(h)` draws a **step line** on the true (possibly variable) bin
edges. It is display geometry without errors and not a fit input.
`render_root(h, GraphSpec(...), "output/hist")` draws it standalone with a visible
zero baseline on linear axes. On log y axes this step path rejects non-positive
contents.

### ROOT histogram draw options (TH1D/TH2D/TH3D)

A standalone `Hist1D` keeps the TGraph step path when `draw` is omitted or `:steps`.
The symbols below are presets; a **string is passed verbatim to `TH1D::Draw`**,
including case and combinations. For example `:points` is the preset `HIST P`,
whereas `"P"` is ROOT's `P`. The empty string requests ROOT's default drawing.
Option semantics follow the [THistPainter documentation](https://root.cern.ch/doc/v640/classTHistPainter.html).

| `draw` | ROOT option | Appearance |
|---|---|---|
| `:hist` | `HIST` | Step outline, no statistical errors |
| `:e`, `:e0`, `:e1` | `E`, `E0`, `E1` | ROOT error bars |
| `:e2` | `E2` | Error rectangles of bin width |
| `:e3` … `:e6` | `E3`–`E6` | Connected error areas and ROOT smoothing variants |
| `:hist_e1` | `HIST E1` | Steps plus error bars |
| `:bar` | `HIST BAR` | Filled bars, no default Sumw2 errors |
| `:points`, `:points_zero` | `HIST P`, `HIST P0` | Bin-center markers |
| `:line`, `:curve`, `:text` | `HIST L`, `HIST C`, `HIST TEXT` | Line, curve, bin-value text |

Errors are the square root of the stored `sumw2`, set with `TH1D::SetBinError`. With
symbol presets, `E2`–`E6` use a translucent fill of the chosen color and `BAR` an
opaque fill at 80% bin width. ROOT's `E4`/`E6` smoothing can exaggerate shapes near
edges; do not read it as a numerical uncertainty. On log y axes, zero contents or
non-positive lower error limits are hidden or clipped by ROOT; bin values are
unchanged.

String options such as `"HBAR3"`, `"PIE"`, `"AXIS"`, `"HIST ]["`, `"E0 E1 X0"`,
`"BAR2"`, `"LF2"`, `"TEXT45"` are forwarded; ROOT then draws the histogram as the
first object. Explicit `GraphSpec` ranges are forwarded as axis ranges/limits;
otherwise ROOT chooses the range. `stats=true` shows ROOT's statistics box; see
[Statistics and fit boxes](#statistics-and-fit-boxes).

`SAME`/`SAMES` are rejected on a standalone canvas and used when drawing several
histograms in order; every later layer must contain `SAME`. `FUNC` requires a valid
`FitResult` passed as `attached_fit`.

```julia
other = Hist1D(h.edges, [1.0, 4.0, 2.0], [1.0, 4.0, 2.0])
fit = fit_root(h, :pol0; initial=[2.0])
render_root(h, spec, "output/pie"; draw="PIE")
render_root([h, other], spec, "output/overlay"; draw=["HIST", "E1 SAME"],
            styles=[GraphStyle(color="#336699"), GraphStyle(color="#cc3311")])
render_root(h, spec, "output/with-attached-fit"; draw="FUNC", attached_fit=fit)
```

`Grid2D` is a native `TH2D`, so `Graph2DSpec(draw="COLZ2")`, `"BOX1"`, `"CANDLEX"`,
or `Graph3DSpec(draw="SURF7Z")` are passed through, and `Graph2DSpec` accepts
`logx`/`logy`/`logz`. 2D overlays use
`render_root([grid1, grid2], Graph2DSpec(contours=8), stem; draw=["COLZ", "CONT3 SAME"])`.
The `:heatmap`, `:contour`, and `:surface` presets remain available.

A true three-axis binned histogram is `Grid3D(xedges, yedges, zedges, values)` with
`values[ix, iy, iz]` copied into `TH3D` bins; it differs from XYZ points (`Graph3D`).

```julia
volume = Grid3D(xedges, yedges, zedges, values)
render_root(volume, Graph3DSpec(draw="BOX2Z"), "output/volume")
render_root(volume, Graph3DSpec(draw="ISO"), "output/isosurface")
```

Option strings select ROOT's allowed combinations; ROOT may silently ignore
options it does not know. Options are limited to 127 bytes of printable ASCII and
may not contain control characters, quotes, backslashes, semicolons, or
unbalanced parentheses. In particular `SPEC` uses Cling inside ROOT, so these
restrictions must not be bypassed.

**Batch PNG/PDF boundary:** OpenGL-only options (`GLLEGO`, `GLSURF`, `GLCOL`,
`GLBOX`, `GLISO`) produce empty batch output and are rejected. `TH2Poly` and
`THStack` are available as `PolygonHist2D` and `HistStack`. Named `TCutG` selection,
interactive picking/projection, and TH3 overlays are not implemented. The fixed
layouts of comparison and fit/pull figures are not changed by `draw`.

```julia
fit = fit_root(h, :pol0; initial=[1.0])
files = render_root(h, fit, GraphSpec(xlabel="Mass", ylabel="Events"),
                    "output/fit-pull"; format=:both)
```

This overload draws **two pads with a shared x axis**: the visible bins and fit
curve on top, and valid χ² pulls below. The default x range is the histogram's
outer edges (override with `GraphSpec.xlimits`). `logy` does not apply to the pull
pad. Specify `pull_limits=(low, high)` or use the automatic range that contains
all valid pulls and zero. `FitResult` does not store the identity of its data, so
pass **the fit obtained from the same `Hist1D`**. Poisson likelihood fits are not
supported in this figure. `scripts/histogram_fit_demo.jl` shows DataFrame
selection → histogram → fit → combined figure.

Arbitrary `TF1` strings are not accepted because they execute code through Cling;
choose from the model list above.

## Statistics and fit boxes

ROOT's statistics box (`TPaveStats`) is drawn by ROOT itself: every number is
computed and formatted by `THistPainter`/`TGraphPainter`, exactly as on a ROOT
canvas. All keywords are opt-in; plots without them are unchanged.

```julia
h = histogram(values, edges).hist
render_root(h, spec, "output/default"; draw="E1", label="signal", stats=true)

fit = fit_root(h, :gaus; initial=[200.0, 0.0, 1.0])
render_root(h, spec, "output/fit"; draw=:e1, stats=true, fitstats=true,
            attached_fit=fit, fit_label="Gaussian fit")

render_root(h, spec, "output/fields"; draw=:hist,
            stats=(:name, :entries, :mean_error, :std_error, :underflow, :overflow),
            stats_position=:topleft)
render_root(xy.hist, Graph2DSpec(), "output/xy"; stats=true)          # Hist2D/Hist3D
render_root(profile, GraphSpec(), "output/profile"; stats=true)       # Profile1D/2D/3D
render_root([graph], GraphSpec(), "output/graph"; attached_fit=gfit, fitstats=true)
PlotPanel(h, spec; draw=:e1, stats=true, attached_fit=fit, fitstats=true)  # in PlotFigure
```

| Keyword | Values | ROOT equivalent |
|---|---|---|
| `stats` | `false` (default), `true` (= name, Entries, Mean, Std Dev), a tuple of fields, a ROOT letter string, or an integer | `gStyle->SetOptStat` |
| `fitstats` | `false` (default), `true` (= χ²/ndf, Prob, parameters ± errors), a tuple of fields, `"pcev"` letters, or an integer | `gStyle->SetOptFit` |
| `attached_fit` | a valid `FitResult` | `TF1` in the object's list of functions |
| `fit_label` | legend label of the attached TF1 (`""` = no legend entry) | |
| `stats_position` | `:auto`, `:topright`, `:topleft`, `:bottomright`, `:bottomleft`, or NDC `(x1, y1, x2, y2)` | `TPaveStats::SetX1NDC`… |
| `stats_name` | box header; default: the label (`"Histogram"`/`"Profile"` when empty) | histogram name |
| `stats_headroom` | `false` (default, ROOT's behaviour) or `true`: raise the automatic y maximum of a 1-D frame so the data stay below a top box | `h->SetMaximum(…)` |

Statistics fields (`stats=(…)`) and their ROOT letters: `:name` (`n`),
`:entries` (`e`), `:mean` (`m`), `:mean_error` (`M`), `:std` (`r`),
`:std_error` (`R`), `:underflow` (`u`), `:overflow` (`o`), `:integral`
(`i`), `:integral_width` (`I`), `:skewness` (`s`), `:skewness_error` (`S`),
`:kurtosis` (`k`), `:kurtosis_error` (`K`). So `stats="nemrou"` equals
`stats=(:name, :entries, :mean, :std, :underflow, :overflow)` and the integer
`110111`. The integer `1` means ROOT's default `1111`, as in ROOT. Fit fields:
`:values` (`v`), `:free_values` (only non-fixed parameters), `:errors` (`e`),
`:chi2` (`c`), `:prob` (`p`); the integer `1` means ROOT's `111`.
`ROOTPlots.stats_option` and `ROOTPlots.fit_option` return the normalized integers.

**What each object shows.**

| Object | Box |
|---|---|
| `Hist1D` (any `draw`, including presets; `:steps` is drawn as `HIST` when a box or fit is requested) | Entries, Mean, Std Dev, …; fit lines with `attached_fit` |
| `Hist2D`, `Hist3D` | Entries, Mean x/y(/z), Std Dev x/y(/z); `:underflow`/`:overflow` print ROOT's flow table |
| `Profile1D` | Entries, Mean, Mean y, Std Dev, Std Dev y |
| `Profile2D`, `Profile3D`, `PolygonHist2D` | ROOT's TProfile2D/TProfile3D/TH2Poly box |
| `GraphSeries` with `attached_fit` | fit box only (`stats` is not accepted; ROOT has no graph statistics) |
| `HistStack` | none: ROOT draws no box for a THStack, so `stats` is rejected |
| Overlays `render_root([h1, h2, …]; draw=[…])` | box of layer 1, plus one per layer drawn with `SAMES` (stacked below, text in the layer color) |

**Entries and moments.** ROOT's Entries is the number of `Fill` calls, and its
Mean/Std Dev use the unbinned sums of the in-range entries. `histogram` records
both (`hist.entries`, including under/overflow fills, and `hist.moments` =
`(Σw, Σw², Σwx, Σwx²[, …])` in `TH1::GetStats` order), plus the flow bins,
which are filled into ROOT's under/overflow bins. The box therefore matches a
ROOT histogram filled event by event. A `Hist1D`/`Hist2D`/`Hist3D` built
directly from bin contents has `entries === nothing` and no moments: ROOT then
shows the **effective entries** `(Σw)²/Σw²` (`TH1::ResetStats`) and bin-centre
moments, as for any ROOT histogram whose fill history is unknown. Pass
`entries=…`, `moments=…`, `underflow=…`, `overflow=…` (or `bins_with_flow=…`
for 2D/3D) to the constructors to supply them. `rebin` keeps them (like
`TH1::Rebin`); `projection` returns plain visible-bin histograms. With
explicit `xlimits` on a raw-option histogram, ROOT recomputes the statistics
from the bins inside the range.

**Fit box.** `attached_fit` attaches the fitted function to the native object
as a `TF1` with the fitted parameters, errors, χ², ndf, and fixed flags, so the
fit box lines come from ROOT (`Prob` is `TMath::Prob(χ², ndf)`). Parameter
names are ROOT's: `Constant`, `Mean`, `Sigma` for `:gaus`; `p0`, `p1`, … for
the other models. For likelihood fits the χ² line shows the stored fit
statistic (`fit.chi2`). The TF1 is drawn by ROOT over `fit.range` with ROOT's
default function style.

**Placement.** `:auto` is ROOT's default (`gStyle` StatX/StatY: top-right
corner of the pad, partly over the frame); on 2D/3D drawings with a color
palette (`Z` options such as `COLZ`) the box moves left of the palette. The
corner symbols place the box inside the frame corner, and an NDC tuple sets it
exactly. A box whose lines would not fit ROOT's default width at ROOT's
line-height text size is widened to the left (keeping its right edge), so long
fields such as `Mean = 0.0282 ± 0.0224` stay readable. An automatic legend
(`legend=:topright` etc., no `legend_spec`) that would overlap the box is moved
below it (or above it for bottom corners); an explicit `legend_spec` is never
moved. In a `PlotFigure` each panel keeps its own options.

**Headroom.** As in ROOT, the default top-right box is drawn over the frame
and may cover data. `stats_headroom=true` applies the usual ROOT remedy
(`SetMaximum`): the y maximum is raised until every bin, drawn error bar,
graph point, and fit-function value in the horizontal span of the box lies
below it; the box position and the y minimum are unchanged, and nothing changes
when the data already clear the box. It applies to 1-D histograms, profiles,
graphs, and 1-D figure panels with a box at `:auto`, `:topright`,
`:topleft` or an NDC box in the upper half of the frame, and it cannot be
combined with explicit `ylimits`.

```julia
render_root(h, spec, "output/fit"; draw=:e1, stats=true, fitstats=true,
            attached_fit=fit, stats_headroom=true)
```

**Testing hook.** If the environment variable `ROOTPLOTS_STATS_DUMP` names a
file, the native renderer appends each box's NDC position, printed lines, and
ROOT's full-precision values (`GetEntries`, `GetMean`, `GetStdDev`, fit
parameters) to it, followed for 1-D frames by a `FRAME` record (y range in pad
coordinates, highest data point under the box and lowest box edge in NDC). The
test suite uses it to compare boxes with Julia values.

## 2D scalar fields and 3D plots

Here **2D** means a value `z(x, y)` shown as color or contours, not the planar
`GraphSeries(x, y)`. `Grid2D`'s `values[ix, iy]` is the value of cell
`xedges[ix:ix+1] × yedges[iy:iy+1]`: the first matrix axis is x and the second is
y, without transposition.

```julia
grid = Grid2D([0., 1., 3.], [0., 2., 5.], [1. 2.; 3. 4.])
render_root(grid, Graph2DSpec(draw=:heatmap, zlabel="Response"), "output/heatmap")
render_root(grid, Graph2DSpec(draw=:contour, contours=8), "output/contours")
render_root(grid, Graph3DSpec(draw=:surface, zlabel="Response"), "output/surface")

cloud = Graph3D([0., 1., 0.5], [0., 0.5, 1.], [1., 2., 1.5])
render_root(cloud, Graph3DSpec(draw=:points), "output/xyz")
```

Grid values are set directly in variable-bin `TH2D` cells. Heatmaps (`COLZ`) show
cell values; contours (`CONT3`) and surfaces (`SURF2Z`) are **renderings between bin
centers**, not exact analytic surfaces or arbitrary meshes. An all-zero grid uses
a display color range of `[-1, 1]` so ROOT does not blank it; values are unchanged.
`Graph3D` keeps duplicate xy points in plain point display. `interpolate=true`
explicitly enables TGraph2D triangulation, which rejects duplicate or collinear xy.

`Graph2DSpec` accepts `draw=:heatmap|:contour` or a ROOT option string,
`contours=1:255`, and `palette=:viridis|:bird`. `Graph3DSpec` accepts `:surface` for
`Grid2D`, `:points` for `Graph3D`, or ROOT option strings for regular histograms, plus
`theta`, `phi`, and `palette`. Both set titles, axis labels, and canvas `size`.
Palettes apply to grids only; XYZ marker color/size come from
`Graph3D(...; style=GraphStyle(...))`. Example: `scripts/graphnd_demo.jl`
(`scripts/graphnd_demo.jl pdf` writes PDFs only).

Grids support linear/log coordinates and TH3 volume options.
`Graph3DSpec(draw=:surface|:triangles|:heatmap|:contour, interpolate=true, resolution=(40, 40))`
is explicit scattered-point interpolation. Log axes for XYZ points are not yet
supported. 3D grid surfaces need at least two cells per axis, and XYZ data need at
least two distinct x and y values. Axis styles (`xaxis`, `yaxis`, `zaxis`; the z
style also formats the palette) and `canvas` are fields of both specs; use
`PlotPanel/PlotFigure` for ranges, color ranges, explicit contour levels, and
multiple pads. 2D/3D fitting
is not implemented; analytic TF2/TF3 drawing is described under
[Analytic TF2/TF3](#analytic-tf2tf3-and-sampled-julia-functions).

## Comparison-panel style and log axes

```julia
spec = PanelSpec("Mass [GeV]"; logx=false, logy=true, ratio_logy=false,
                 ymin=0.1, ymax=10_000.0, ratio_limits=(0.5, 1.5))
style = PlotStyle(marker_style=20, marker_size=0.9, marker_color="#000000", line_width=2,
                  stack_fill_style=1001, band_fill_style=3345, band_color="#000000",
                  legend_columns=2)
render_root([cmp], [spec], "output/mass"; style)
```

`logx` maps to `SetLogx` on both pads, `logy` to the upper pad, and `ratio_logy` to
the lower pad. `ymin`/`ymax` and `ratio_limits` are **physical coordinates, not
logarithms**; when omitted, log ranges are derived from drawable positive values.
Log x requires all bin edges to be positive. Zero points are hidden on log y,
while error bars and bands extending into the positive range are clipped. The
numerical `Hist1D`/`Comparison` data are not modified.

`PlotStyle` applies to the whole comparison figure; codes go directly to ROOT's
`SetMarkerStyle`/`SetFillStyle`. Unspecified values keep the renderer defaults.
Component colors are set per `Component`. When a panel omits its ratio range, the
`ratio_limits` keyword of `render_root` is used. `PanelSpec`'s `annotation` and
`lumi` are free text drawn above the frame; no experiment label is added.

## Experiment styles

`style` selects a ROOT `TStyle` for the canvas, and `cms_label` / `atlas_label` draw
an experiment label with ROOT `TLatex`. Both are opt-in: without them nothing
changes.

| Option | Specs | ROOT |
|---|---|---|
| `style=:Name` | `GraphSpec`, `Graph2DSpec`, `Graph3DSpec`, `PanelSpec`, `PlotFigure` | `gROOT->SetStyle("Name")` before any canvas or histogram exists |
| `cms_label` | the same specs except `PlotFigure` | `TLatex`: "CMS" (font 61), extra text (font 52), luminosity/energy (font 42) |
| `atlas_label` | the same specs except `PlotFigure` | `TLatex`: "ATLAS" (font 72), text and lines (font 42) |

**Styles.** Names are matched case-insensitively against ROOT's style list
(`gROOT->GetListOfStyles()`); `root_styles()` returns them, e.g. `"ATLAS"`,
`"BELLE2"`, `"Bold"`, `"Classic"`, `"Default"`, `"Modern"` (ROOT's default), `"Plain"`,
`"Pub"`, `"Video"`, plus `"CMS"`. An unknown name fails in the native backend with the
list of available names. No `gROOT->ForceStyle()` is used: nothing is read from a
file, and every object is created after the style is selected, so it takes the
style's attributes, while explicit settings applied afterwards keep their values.

- Explicit options win: `AxisStyle`, `CanvasStyle` (margins, ticks, grids, title box),
  `LegendSpec`, statistics/fit-box options, and draw options are applied after the
  style. The package's opt-in boxes stay off (`SetOptStat(0)`, `SetOptFit(0)`,
  `SetOptTitle(0)`) until requested, whatever the style says.
- Package layout defaults yield to the style where they are cosmetic: flat 2-D
  views (Hist2D, Grid2D, profiles, TF2) keep the style's left/bottom/top margins
  (a palette still gets room on the right: at least 0.18, or 0.2 with a title, widened
  until the palette labels, `×10^n` exponent and z title fit at the style's z label and
  title sizes, with a z title offset clear of the labels; a top margin too small for
  the exponent grows as well);
  automatic legends use 0.8 x the style's label size; the fit/pull and
  comparison views take the style's side margins, ticks, and text sizes,
  rescaled to each pad's height so text keeps its one-pad size (see "Two-pad
  views" below), and the comparison legend moves into the frame. 3-D views keep
  the package's margins.
- The style is process-global (one per canvas): figure panels and comparison
  panels that name a style must agree. `PlotFigure(...; style=:CMS)` gives the
  style to every panel spec without one.

**CMS.** ROOT has no built-in CMS style. `style=:CMS` creates a `"CMS"` `TStyle` in the
backend from ROOT setters, with values that follow the public CMS publication
style conventions: white canvas, pads and frame without borders; no statistics,
title, date or file boxes; Helvetica (font 42) for all text; axis title size
0.06 (x offset 0.9, y offset 1.25), label size 0.05 (offset 0.007); ticks on all
four frame sides (`SetPadTickX/Y(1)`), tick length 0.03, 510 divisions; pad margins
top 0.08, bottom 0.13, left 0.16, right 0.04; markers 20, function line width 2,
end-error size 2. The palette is the package's `palette` option (viridis by
default). In the comparison view, y-axis numbers use an exponent above 999
(`TAxis::SetMaxDigits(3)`; `yaxis=AxisStyle(max_digits=...)` overrides).

The single switch is `style=:CMS`: it also draws the default label
"CMS Preliminary". `cms_label` refines or removes it:

| `cms_label` | Effect |
|---|---|
| `nothing` (default) | label "CMS" + "Preliminary" when `style=:CMS`, otherwise none |
| `false` / `true` | no label / default label (with any style) |
| `(extra=..., lumi=..., energy=..., position=...)` | fields below |

- `extra`: text after "CMS" (default `"Preliminary"`; e.g. `"Simulation"`, `"Private Work"`;
  `""` for none), in italics (font 52), 0.76 x the "CMS" size.
- `lumi`, `energy`: text above the frame on the right, `"lumi (energy)"`, e.g.
  `cms_label=(lumi="138 fb^{-1}", energy="13 TeV")` gives "138 fb^{-1} (13 TeV)";
  `energy` alone gives "(13 TeV)". A comparison panel's `lumi` is used when neither
  is given.
- `position`: `:inside_left` (default for 1-D plots; "CMS" bold, font 61, inside the
  frame at the top left with the extra text below it), `:inside_right`, or `:outside`
  ("CMS" and the extra text on one line above the frame at the left; the default
  for 2-D/3-D views so the label does not cover coloured cells).
- Sizes follow the CMS convention of scaling with the pad's top margin: "CMS" is
  0.75 x the top margin, the luminosity text 0.6 x. In figure pads they scale with
  each pad; in the upper pad of a two-pad view they keep their one-pad size.

**ATLAS.** `style=:ATLAS` selects ROOT's built-in `"ATLAS"` style unchanged.
`atlas_label=(text="Internal", energy="13 TeV", lumi="140 fb^{-1}")` draws
"ATLAS" in bold italics (font 72) followed by the text (font 42, default
`"Internal"`; e.g. `"Preliminary"`), and a second line "#sqrt{s} = 13 TeV, 140 fb^{-1}"
below it, at the style's text size (`gStyle->GetTextSize()`, 0.05 in the ATLAS
style). `position=:inside_left` (default), `:inside_right`, or `:outside` (the
second line goes above the frame on the right; when the row is too narrow for
both, it is stacked under "ATLAS <text>" at the left). The ATLAS label is opt-in;
`atlas_label=true` draws "ATLAS Internal".

**Room above the frame.** Label texts drawn above the frame (outside labels, the
CMS luminosity text) always end inside their pad: unless `CanvasStyle(top_margin=...)`
fixes it, the pad's top margin grows until they fit with a small white edge
(e.g. the outside ATLAS label under the ATLAS style's 0.05 top margin, in single
canvases, figure pads and the upper pad of two-pad views). With a fixed top margin
too small for them, the texts move down instead.

**Comparison panels** (`PanelSpec`): with a style or a label, both pads span the full
canvas height and take the style's margins and text sizes as described above. A
label carries the panel texts: `lumi` is the CMS luminosity text (or the ATLAS
second line) and `annotation` a line below the label. Without a label, under a
style, `annotation` and `lumi` are drawn above the frame (left and right).

**Fit/pull view**: the label goes to the upper pad.

The CMS and ATLAS labels are for collaboration use under each collaboration's
publication rules; the package only provides the styling.

Testing hook: `ROOTPLOTS_AXIS_DUMP` also records the selected `gStyle` (`GSTYLE` row:
name, pad margins, label/title font and size, pad ticks, `OptStat`) and every
`TLatex` in each pad (`LATEX` rows: text, font, size, position, alignment, NDC flag);
see `test/experiment_style_tests.jl` and `scripts/experiment_style_demo.jl`.

## Backend architecture

```text
vectors/DataFrame → Julia numeric objects → versioned TSV → isolated C++ ROOT process → PNG/PDF and fit results
```

The two native bridges are compiled with ROOT's `root-config --cflags --libs` and
stored in the Julia depot's Scratch cache; the package source is never written.
Paths and values are passed as argument arrays and an explicit text protocol,
never evaluated as ROOT macros. No `.root` files are created. ROOT's `TLatex`
interprets `#`, `_`, and `^`; `μ` and `ℓ` are converted to `#mu` and `#ell`. Pixel
identity across Unicode text, fonts, or operating systems is not guaranteed.

**ROOT messages.** Each native process sets ROOT's `gErrorIgnoreLevel` to
`kWarning`, so routine `Info in <TCanvas::Print>: ... has been created` lines
are not printed, while ROOT warnings and errors still reach stderr. With the
environment variable `ROOTPLOTS_ROOT_VERBOSE=1` the processes keep ROOT's
default message level (e.g. `withenv("ROOTPLOTS_ROOT_VERBOSE" => "1") do ... end`).

## Axes, titles and number format

ROOTPlots is a Julia front end to ROOT's own drawing, so every option below is
one ROOT setter with a recognisable name, applied by ROOT to the drawn axis or
pad. Leaving a field at `nothing` keeps ROOT's default (or the documented layout
default of the painter); nothing is reformatted in Julia. All options are
opt-in: existing plots are unchanged.

```julia
spec = GraphSpec(xlabel="Energy (keV)", ylabel="Events", title="Spectrum", logx=true, logy=true,
    xaxis=AxisStyle(more_log_labels=true, title_offset=1.3),
    yaxis=AxisStyle(no_exponent=true, label_color="#cc3311"),
    canvas=CanvasStyle(grid_x=true, grid_y=true, opt_title=true, title_x=0.15, title_align=13))
render_root(h, spec, "output/spectrum"; draw="E1")
```

### AxisStyle and LabelChange

| field | ROOT call | notes |
|---|---|---|
| `title_font`, `title_size`, `title_offset`, `title_color` | `TAxis::SetTitleFont/Size/Offset/Color` | |
| `label_font`, `label_size`, `label_offset`, `label_color` | `TAxis::SetLabelFont/Size/Offset/Color` | |
| `axis_color` | `TAxis::SetAxisColor` | axis line and ticks |
| `tick_length` | `TAxis::SetTickLength` | |
| `ticks` | `TAxis::SetTicks("+" / "-" / "+-")` | tick side; `CanvasStyle.ticks_x/y` add ticks on the opposite frame edge (`TPad::SetTickx/y`) |
| `divisions`, `optimize_divisions` | `TAxis::SetNdivisions(n, optim)` | `optimize_divisions=false` keeps exactly `n` (ROOT stores `-n`) |
| `max_digits` | `TAxis::SetMaxDigits` | |
| `no_exponent` | `TAxis::SetNoExponent` | 400000 instead of 400 ×10³ |
| `more_log_labels` | `TAxis::SetMoreLogLabels` | label 2, 3, 5, … on short log ranges |
| `decimals` | `TAxis::SetDecimals` | 1.0 instead of 1 next to 0.5 |
| `center_title`, `center_labels` | `TAxis::CenterTitle`, `TAxis::CenterLabels` | |
| `time_display`, `time_format` | `TAxis::SetTimeDisplay`, `TAxis::SetTimeFormat` | strftime format; ROOT's `%F` offset syntax is accepted |
| `time_offset` | `TAxis::SetTimeOffset(t, "gmt")` | UNIX time (s) of axis value 0, always GMT for reproducible output |
| `exponent_offset` | `TGaxis::SetExponentOffset(dx, dy, "x"/"y")` | static in ROOT: one value per direction per canvas (conflicts are rejected) |
| `bin_labels` | `TAxis::SetBinLabel(bin, label)` | histogram bin axes only; strings (bins 1, 2, …), `bin => label` pairs, or a Dict |
| `labels_option` | `TAxis::LabelsOption("h"/"v"/"u"/"d")` | orientation of bin labels; needs `bin_labels` |
| `change_labels` | `TAxis::ChangeLabel` | vector of `LabelChange(index; angle, size, align, color, font, text)` |

`LabelChange(index; ...)` follows ROOT's `ChangeLabel(labNum, labAngle, labSize,
labAlign, labColor, labFont, labText)`: `index` 1 is the first drawn label and -1
the last; unset fields keep the label's attribute; `size=0` hides a label.
`ChangeLabel` acts on numeric labels, not on bin labels, as in ROOT.

Bin labels name histogram bins, so they are accepted where the drawn axis is a
histogram axis: Hist1D (all draw modes; a default `:steps` plot then draws the
native TH1D as `HIST`, and an owned frame takes the histogram's binning), Hist1D
overlays and Hist1D+graph plots, stacks, profiles, Hist2D/Hist3D/Grid2D/Grid3D
axes, the first Hist1D layer of a 1-D panel, and the comparison x axis. Graphs,
functions, XYZ graphs, TF2/TF3, TH2Poly and special graphs have no bin axis and
reject them (use `change_labels` there). Bin indices are checked against the
bin count before ROOT runs.

### CanvasStyle: grids and the ROOT title box

| field | ROOT call |
|---|---|
| `grid_x`, `grid_y` | `TPad::SetGridx/SetGridy` |
| `opt_title` | `gStyle->SetOptTitle`: ROOT draws the plot title as its `TPaveText` "title" |
| `title_x`, `title_y`, `title_w`, `title_h` | `gStyle->SetTitleX/Y/W/H` (NDC; 0 = ROOT's automatic size) |
| `title_align` | `gStyle->SetTitleAlign` (ROOT default 23: centred at the top) |
| `title_font`, `title_font_size` | `gStyle->SetTitleFont(font, "t")`, `gStyle->SetTitleFontSize` |
| `title_border_size`, `title_fill_color`, `title_style`, `title_text_color` | `gStyle->SetTitleBorderSize/SetTitleFillColor/SetTitleStyle/SetTitleTextColor` |

By default (`opt_title` unset) a non-empty title is the package's TLatex heading
at the top left, exactly as before. With `opt_title=true` the title string goes to
`TH1::SetTitle` of the frame, so ROOT's `"title;x title;y title"` form also sets the
axis titles. The `title_*` fields require `opt_title=true`. gStyle is process-global:
the values are set only while ROOT creates this pad's title box, and the box is
then kept (`TH1::kNoTitle` on the frame), so each pad of a figure has its own
title settings. Polar and pie charts have no Cartesian frame and reject axis
styles, grids, and the title box.

### 2D/3D, panels, and two-pad views

- `Graph2DSpec`/`Graph3DSpec` have `xaxis`, `yaxis`, `zaxis`, and `canvas`. The z style of
  a 2-D view also formats the palette, because ROOT's palette axis imports the
  histogram's z-axis attributes. Explicit `canvas` margins override the painter's
  layout margins (e.g. the palette margin); an explicit z `title_offset` is kept as given.
- `PlotPanel` keeps its `xaxis`/`yaxis`/`zaxis`/`canvas` keywords; for 2D/3D panels the
  spec's fields are the defaults the keywords override. Nested figures apply the
  same settings per pad.
- `PanelSpec` (comparison with ratio pad): `xaxis` styles the shared x axis (drawn
  on the ratio pad; the upper pad keeps hidden labels but the same divisions/time
  settings), `yaxis` the upper y axis, `ratio_yaxis` the ratio y axis, and
  `canvas` accepts `grid_x`, `grid_y`, `ticks_x`, `ticks_y` for both pads.
- Fit/pull figure (`render_root(h, fit, spec, stem)`): `xaxis` applies to both pads,
  `yaxis` to the upper pad, grids/ticks to both pads, and the title box to the
  upper pad. Its margins, annotations, and legend stay layout-owned.

**Two-pad views: text sizes and edge labels.** ROOT gives axis text sizes as a
fraction of the smaller pad dimension (pixels for precision-3 fonts, e.g. 43), so
one value is larger in the short ratio/pull pad than in the upper pad. The
comparison and fit/pull layouts therefore size the lower pad from the upper y
axis, after the request's axis styles:

- The lower pad's y labels and title and its x labels and title get the upper y
  axis's label and title sizes in pixels (under a style: the style's one-pad
  sizes, e.g. 0.05/0.06 of the canvas width under `style=:CMS`). A lower y
  title longer than its pad (e.g. "Data / pred." at the CMS size) is reduced
  until it fits.
- The lower y title is placed at the upper title's distance from the axis, so
  the two titles line up (`TAxis::SetTitleOffset` is converted between the pads,
  since TGaxis places a y title 1.6 x size x offset pad widths from the axis). ROOT's
  automatic y title offset (0, e.g. ROOT's default style) is computed per pad from
  its own labels; in a two-pad view the upper title is then placed next to the
  left pad edge instead.
- The upper pad's bottom margin and the lower pad's top margin leave room for
  half a y label, so labels at the frame edges (1.5 and 0.5 of a 0.5-1.5 ratio
  axis, 0 of a linear upper axis) are drawn whole inside their pads. To drop an
  edge label instead, use `ratio_yaxis=AxisStyle(change_labels=LabelChange(-1; size=0))`
  or `divisions`.
- Explicit `AxisStyle` sizes and offsets win and keep ROOT's meaning for the pad
  they are drawn in: `yaxis` sizes are fractions of the upper pad, `ratio_yaxis`
  and `xaxis` sizes fractions of the ratio pad. Sizes left unset in the lower
  pad still follow the upper y axis, so `yaxis=AxisStyle(title_size=0.035,
  label_size=0.035)` makes the whole figure's axis text that size. For one
  pixel size given directly, use a precision-3 font, e.g.
  `AxisStyle(title_font=43, title_size=24)` (24 pixels in either pad).

### Coverage (audit)

The table records, per entry point, the state before this feature → now.
"AxisStyle" means the full field set above; "basic" is the earlier subset
(fonts, sizes, offsets, tick length, divisions, max digits, centred title).

| entry point | log x / y / z | plot title | axis titles | AxisStyle | CanvasStyle (grids, title box) |
|---|---|---|---|---|---|
| graphs `render_root(series, GraphSpec)` | x, y | TLatex → TLatex or ROOT box | x, y | basic x/y → x/y | no grids → yes |
| Hist1D presets, raw options, `:steps` | x, y | TLatex → TLatex or ROOT box | x, y | basic → full (bin labels included) | → yes |
| Hist1D overlays, Hist1D + graphs | x, y | TLatex → TLatex or ROOT box | x, y | basic → full | → yes |
| mixed layers / bent & multi-error graphs (panel) | x, y | TLatex → TLatex or ROOT box | x, y | basic → full | → yes |
| `ModelCurve` (TF1) | x, y | TLatex → TLatex or ROOT box | x, y | basic → full (no bin labels) | → yes |
| `HistStack` | x, y | TLatex → TLatex or ROOT box | x, y | basic → full | → yes |
| `Profile1D` | x, y | TLatex → TLatex or ROOT box | x, y | basic → full | → yes |
| `Profile2D`/`Profile3D` | x, y, z | TLatex → TLatex or ROOT box | x, y, z | none → x/y/z | none → yes |
| `PolygonHist2D` (TH2Poly) | x, y, z | TLatex → TLatex or ROOT box | x, y, z | none → x/y/z (no bin labels) | none → yes |
| `Grid2D` heatmap/contour/raw, overlays | x, y, z | TLatex → TLatex or ROOT box | x, y, z | none → x/y/z (+palette) | none → yes |
| `Grid2D` surface (`Graph3DSpec`), `Grid3D` (TH3) | x, y, z | TLatex → TLatex or ROOT box | x, y, z | none → x/y/z | none → yes |
| `Hist2D`/`Hist3D` (one panel) | x, y, z | TLatex → TLatex or ROOT box | x, y, z | basic via panel keywords → full, also from the spec | panel keyword → also from the spec |
| `Graph3D`/`Graph3DErrors` (XYZ) | rejected (ROOT's TGraph2D painter) | TLatex → TLatex or ROOT box | x, y, z | none → x/y/z (no bin labels) | none → yes |
| `ModelSurface` (TF2) | x, y, z | TLatex → TLatex or ROOT box | x, y, z | none → x/y/z (no bin labels) | none → yes |
| `ModelVolume` (TF3) | rejected (batch TF3 painters) | TLatex → TLatex or ROOT box | x, y, z | none → x/y/z (no bin labels) | none → yes |
| QQ, spline, efficiency | x, y | TLatex → TLatex or ROOT box | x, y | basic → full (no bin labels) | → yes |
| polar, pie | not applicable | TLatex | n/a | silently ignored → rejected | grids/title box rejected |
| `PlotPanel`/`PlotFigure`, nested | per panel as above | per pad | per pad | per pad | per pad |
| comparison (`PanelSpec`) | x (both pads), y, ratio y | annotation/lumi text | x, y | none → x, y, ratio y | none → grids, ticks |
| fit/pull two-pad | x (both), y (upper) | TLatex → TLatex or ROOT box | x, y | rejected → x (both pads), y | rejected → grids, ticks, title box |

Log scales were already available wherever ROOT's batch painters support them;
the remaining exceptions are ROOT limitations (XYZ graph points/errors and TF3),
and polar/pie charts, whose axes are not Cartesian.

### Not exposed, and why

- **Reversed axes**: ROOT 6.40 has no reversed-axis option for TH1 painting
  (`"RX"` is ignored by THistPainter); TGraph's `"RX"/"RY"` work only when the graph
  draws its own axes, while ROOTPlots draws graphs on an owned frame.
- **`SetRangeUser`**: covered by `xlimits`/`ylimits` (GraphSpec), `ymin`/`ymax`/`ratio_limits`
  (PanelSpec), and the `xlimits`/`ylimits`/`zlimits` keywords of 2D/3D panels
  (`render_root(h2, spec, stem; xlimits=...)` forwards them).
- **Global `TGaxis::SetMaxDigits`**: the per-axis `TAxis::SetMaxDigits` (`max_digits`) does the
  same without leaking into other pads.
- **Sorting `LabelsOption` ("a", ">", "<")**: reorders bin contents, a data
  operation rather than a display option. `ChangeLabelByValue` is not exposed.
- **`logz` in CanvasStyle**: z scales belong to `Graph2DSpec`/`Graph3DSpec` (`logz`);
  1-D views have no z axis.
- The pull axis of the fit/pull figure and the comparison layout's margins and
  title area are layout-owned.

### Testing hook

With the environment variable `ROOTPLOTS_AXIS_DUMP=<file>` both native bridges
append, just before the image is saved, each pad's ROOT state (`GetLogx/y/z`,
`GetGridx/y`, `GetTickx/y`, margins, the "title" box text and attributes) and the
frame axes as ROOT reports them (`GetNdivisions`, `GetMaxDigits`, `GetNoExponent`,
`GetMoreLogLabels`, `GetDecimals`, `GetTicks`, `GetTimeDisplay`/`GetTimeFormat`, colours,
sizes, bin labels, `ChangeLabel` entries). `test/axis_format_tests.jl` checks these
values; `scripts/axis_format_demo.jl` draws every option.
