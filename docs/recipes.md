# Recipes

[Quick start](../README.md) · [API](api.md) · [Testing](testing.md)

Task-based examples. All data are synthetic. Each recipe states what the result
means statistically; the [API](api.md) has the full rules.

## Combine independent samples with per-sample weights

Scale each sample, then merge. Scaling multiplies every bin (flows included) by
the factor and every variance by its square; merging adds bins, variances,
entries, and unbinned moments. The result equals one histogram filled with the
per-event weights, so ROOT's statistics box shows the same Entries, Mean,
Std Dev, Underflow and Overflow as a direct fill.

```julia
using ROOTPlots
edges = collect(-4.0:0.2:4.0)
signal = histogram(x_signal, edges).hist          # generated events, weight 1
background = histogram(x_background, edges).hist
lumi = 140.0                                      # fb^-1
total = merge((lumi * 0.12 / 50_000) * signal,    # cross section [fb] / generated events
              (lumi * 3.4 / 200_000) * background)
render_root(total, GraphSpec(xlabel="m [GeV]", ylabel="Events / 0.2"), "output/total";
            draw=:hist, stats=(:entries, :mean, :std, :underflow, :overflow))
```

- Only **independent** samples may be merged. `merge(h, h)` is rejected; to count
  one sample twice use `2 * h` (variance × 4, not × 2). Copies of overlapping
  samples are not detectable: keep track of which events each histogram holds.
- Different binnings are rejected; `rebin(h, edges)` both histograms onto common
  edges first.
- Entries stay fill counts (scaling does not change them). If one input has
  unknown entries or moments (a histogram built from bin contents), the merged
  entries or moments are unknown and ROOT shows effective entries or bin-centre
  moments instead of mixing the two kinds.
- Signed weights add as they are. Bins can cancel while `Σw²` keeps growing; the
  effective entries `(Σw)²/Σw²` then describe the signed sum, not a number of
  events.
- There is no normalization or subtraction helper. Dividing by a measured
  integral correlates the bins, and a difference of correlated samples does not
  have variance `Σw²₁ + Σw²₂`; do such steps explicitly with the correct
  uncertainty model.

`scripts/histogram_ops_demo.jl` draws the merged and the directly filled histogram
side by side; their ROOT statistics boxes are identical.

## Project a 2-D or 3-D histogram

```julia
h2 = histogram((x, y), (xedges, yedges)).hist
visible = projection(h2, :x)               # visible bins only, statistics unknown
allx = projection(h2, :x; flow=true)       # every event; y flows included, x flows kept
xz = projection(h3, (:x, :z); flow=true)   # x on the x axis, z on the y axis
```

| Want | Use |
| --- | --- |
| The x distribution of all events (like a direct fill of x) | `flow=true` |
| Only events inside the y range | the default (`flow=false`), or restrict before filling |
| ROOT's `ProjectionX()` default range | `flow=true` (contents agree; statistics may differ, see below) |

With `flow=true` the original entries and the matching unbinned moments are kept
only when they are exact: no integrated-axis flow cell under a visible result bin
may hold any weight (neither `Σw` nor `Σw²`). Otherwise the statistics are
unknown and ROOT shows bin-centre Mean/Std Dev with effective entries. The rule is
stricter than ROOT for cancelling signed weights, and keeps exact statistics
where ROOT 6.40 falls back to bin centres because the projection has flow bins.
The default `flow=false` is unchanged from earlier releases and never carries
statistics.

`flow=true` relies on stored flow bins. Histograms from `histogram` have them; a
`Hist2D`/`Hist3D` built from visible arrays without `bins_with_flow` is treated
as having empty flow bins, so its stored entries and moments are reused even if
the original events had flows. Pass `bins_with_flow` (or leave `entries` and
`moments` unset) when the flows are not known to be empty.

## Read fit errors: weighted and unweighted fits

```julia
x, y = collect(0.0:4.0), [1.1, 2.9, 5.2, 6.8, 9.1]
weighted = fit_root(GraphSeries(x, y; yerrors=fill(0.2, 5)), :pol1; initial=[0.0, 1.0])
plain = fit_root(GraphSeries(x, y), :pol1; initial=[0.0, 1.0])
```

| | `weighted` (errors given) | `plain` (no errors) |
| --- | --- | --- |
| parameters | 1.04, 1.99 | 1.04, 1.99 |
| covariance | `inv(X'WX)` = [0.024 −0.008; −0.008 0.004] | `(RSS/ndf)·inv(X'X)`, RSS = 0.107, ndf = 3 |
| `chi2` | 2.675, a goodness-of-fit statistic | 0.107 = RSS with unit errors, not a goodness of fit |
| `normalized_errors` | `false` | `true` |

- With known errors the parameter errors come from the χ² curvature alone; they
  do not shrink or grow with the fit quality. Check `chi2/ndf` and
  `chi2_pvalue(fit)`.
- Without errors ROOT fits with unit errors and scales the covariance by
  `χ²/ndf` (`normalized_errors=true`): the errors then estimate the scatter from
  the residuals, and `chi2_pvalue` refuses such a fit.
- Histogram fits use bin centres and `sqrt(Σw²)` errors; bins with zero error
  (empty bins in particular) are left out of a χ² fit. For low counts use
  `objective=:poisson` (unweighted counts) or `:weighted_poisson`.
- Graph x errors: `xerrors=:effective` (default) uses ROOT's effective variance
  `σ_y² + (f′(x) σ_x)²`; `xerrors=:ignore` (ROOT option `EX0`) fits with y errors only.
- `fit.status == 0 && fit.valid` means the minimizer converged, not that the model
  is right.

## CMS-style figure with a long y-axis title

The `style=:CMS` default follows cmsstyle's `cmsCanvas` geometry; on the default
800×600 canvas its y title offset 0.78 suits short titles only, and a title such
as "Events / GeV" reaches the left canvas edge. Use cmsstyle's own recipe for
longer titles, `cmsCanvas(extraSpace=0.05, yTitOffset=1.1)`:

```julia
spec = GraphSpec(xlabel="m_{ll} [GeV]", ylabel="Events / GeV", style=:CMS,
                 cms_label=(lumi="138 fb^{-1}", energy="13 TeV", extra_space=0.05, y_title_offset=1.1))
```

## Reproduce a figure later

```julia
files = render_root(hist, fit, GraphSpec(xlabel="m [GeV]"), "out/fit"; provenance=true)
# files.request == "out/fit.request.tsv", files.provenance == "out/fit.provenance.toml"
```

Keep both files with the figure. `out/fit.request.tsv` contains the plotted
numbers and labels verbatim; share it only where the data may go. To redraw (same
package version and ROOT):

```julia
using ROOTPlots, TOML
record = TOML.parsefile("out/fit.provenance.toml")
built = build_root()
exe = record["backend"]["executable"] == "rootplots-render" ? built.renderer : built.graphfit
run(`$exe out/fit.request.tsv out/fit-replay --formats=png,pdf`)
bytes2hex(open(ROOTPlots.sha256, "out/fit-replay.png")) == record["outputs"]["png"]["sha256"]
```

## Bound a batch of ROOT jobs

```julia
results = try
    rootmap(samples; workers=4, timeout=300) do s
        render_root(s.hist, GraphSpec(), "out/" * s.name; format=:pdf)
    end
catch e
    e isa RootMapError || rethrow()
    for (i, failure) in e.failures
        failure.ex isa RootBackendError && failure.ex.timedout && @warn "sample $i timed out"
    end
    rethrow()
end
```

A job whose ROOT process takes longer than 300 s fails alone; Ctrl-C stops the
whole batch and its ROOT processes, and leaves other work alone.
