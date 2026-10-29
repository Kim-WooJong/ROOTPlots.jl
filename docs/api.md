# ROOTPlots API and computational semantics

[Quick start](../README.md) · [Installation](installation.md) · [Migration](migration.md) · [Testing](testing.md)

This document describes the supported scope and the scientific contracts of the
computations. ROOTPlots is not a complete ROOT API binding.

## Supported plots

| Task | Julia API | Native ROOT object / behavior |
|---|---|---|
| 1D histograms: errors, bars, steps, overlays | `Hist1D`, `histogram`, `render_root` | TH1D, ROOT draw options |
| Weighted 2D/3D histograms with errors | `Hist2D`, `Hist3D` | TH2D/TH3D contents and Sumw2 errors |
| Heatmaps, contours, LEGO, surfaces, volumes | `Grid2D/3D`, `Graph2DSpec/Graph3DSpec` | TH2D/TH3D painters |
| Stacked / unstacked histograms | `HistStack` | THStack; data/model comparison panels |
| 1D/2D/3D profiles | `Profile1D/2D/3D`, `profile_stats` | TProfile/TProfile2D/TProfile3D filled natively |
| Polygon bins | `PolygonHist2D` | TH2Poly (convex polygons without holes) |
| XY points, lines, curves, fills, bars | `GraphSeries`, `GraphStyle` | TGraph and native draw options |
| Symmetric/asymmetric errors, bands, boxes, arrows | `GraphSeries(...; xerrors, yerrors, errors)` | TGraphErrors/TGraphAsymmErrors |
| XYZ points/lines, errors, explicit interpolation | `Graph3D`, `Graph3DErrors` | TGraph2D/TGraph2DErrors |
| Sampled Julia functions / native TF1 | `function_curve` / `ModelCurve` | Sampled TGraph / TF1 of a predefined model |
| Polar, QQ, spline | `PolarGraph`, `QQGraph`, `SplineCurve` | TGraphPolar, TGraphQQ, TSpline3/5 |
| Pie charts, binomial efficiency | `PieChart`, `EfficiencyGraph` | TPie, TEfficiency with native intervals |
| Lines, LaTeX text, legends, axes, margins | `PlotLine`, `PlotText`, `LegendSpec`, `AxisStyle`, `CanvasStyle` | TLine/TLatex/TLegend/ROOT axes |
| Arrows, boxes, ellipses, secondary axes | `PlotArrow`, `PlotBox`, `PlotEllipse`, `SecondaryAxis` | TArrow/TBox/TEllipse/TGaxis |
| Mixed multi-panel figures, vector PDF | `PlotPanel`, `PlotFigure` | One TCanvas with TPads (not stitched PNGs) |
| Parallel work | `rootmap`; `histogram(...; threads=...)` | Julia tasks/CPU threads + isolated ROOT processes |

Out of scope or not implemented: ROOT file I/O, RDataFrame, GUI/OpenGL,
RooFit/RooStats, arbitrary C++ or TF1 expression strings, arbitrary analytic
TF2/TF3, TGraphBentErrors/TGraphMultiErrors, non-convex or holed TH2Poly bins,
arbitrary nested pads, and geometry/event displays. Not every axis detail of
polar/XYZ plots and not every ROOT class method is exposed. A Julia function can
be sampled and overlaid, but that is not described as the corresponding special
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
displays pre-aggregated values for convex polygons and does not fill events.

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
panels may mix `Hist1D`, `GraphSeries`, and `ModelCurve` layers; ND panels
draw one ND object. The `PlotPanel` keywords `xaxis/yaxis/zaxis`,
`xlimits/ylimits/zlimits`, `canvas`, `lines/texts/legend_spec`,
`style/styles`, and `label/labels` apply to basic XY/ND panels. Advanced,
ModelCurve, XYZ, and special objects use their own spec and style; common
`primitives` may be added to any panel. Because a basic XY panel owns its axes,
use standalone `render_root` for special TH1 options such as HBAR/PIE/FUNC/AXIS.

Palettes and error-cap/PostScript line scale are ROOT global settings, so
conflicting values within one figure are rejected; ND panels must share a
palette. Nested layouts and automatic shared axes are not provided. For grid
panels, `zlimits` is the color range for TH2 and the spatial z range for TH3.
Arrow/Box/Ellipse/SecondaryAxis use data coordinates; a secondary axis is a
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

Runnable examples: `scripts/graphics_gallery.jl` (9 panels),
`scripts/special_graph_demo.jl` (6 panels), `scripts/advanced_histogram_demo.jl`,
`scripts/graph_extension_demo.jl`, `scripts/edge_profile_demo.jl`. All inputs
are synthetic.

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
  `center_title`. Omitted fields keep ROOT defaults; precision-3 fonts (e.g. 43)
  use pixel sizes.
- `CanvasStyle`: four margins, ticks on both sides, `frame_line_width`,
  `end_error_size`, `line_scale_ps`. Margins are canvas fractions.
- `PlotLine` uses **data coordinates** by default; `ndc=true` uses 0–1 canvas
  coordinates. A line `label` is collected into the automatic legend.
- `PlotText` uses **NDC coordinates** by default; `ndc=false` uses data
  coordinates. `font`, `size`, `color`, `align`, and `angle` are configurable.
  Math uses ROOT TLatex syntax (`#sigma`, `#pm`, `#mum`, `mm^{2}`), not a full
  LaTeX engine.
- `legend_spec` replaces the `legend=:topright`-style position. `background=nothing`
  is transparent; a color string is opaque. `entries=nothing` collects labels
  automatically and `entries=[]` shows none. `LegendEntry` sets order, target, and
  the `l/p/pe/ep/lp/lep/f` marker.
- Annotations never extend the data range. Data-coordinate lines/text must be
  positive on log axes. Clipping of out-of-range annotations follows ROOT pads.
- These options apply to **single-pad GraphSeries/Hist1D/TH1 overlays**. Comparison
  panels use `PanelSpec/PlotStyle`. The two-pad fit/pull figure rejects unsupported
  scene settings instead of silently ignoring them; 2D/3D plots use their own specs.

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
moments cannot be recovered from a `Hist1D`. `histstats` bundles these values and
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
otherwise ROOT chooses the range. `stats=true` shows a bin-center mean/std box
but **no Entries**, since the original fill count is unknown.

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
least two distinct x and y values. Use `PlotPanel/PlotFigure` for ND axis styles,
ranges, color ranges, explicit contour levels, and multiple pads. 2D/3D fitting
and arbitrary analytic TF2/TF3 are not implemented.

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
