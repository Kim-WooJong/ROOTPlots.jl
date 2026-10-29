# ROOTPlots.jl

**Julian plotting and fitting with native CERN ROOT.**

Create scientific plots and perform ROOT fits from Julia. Numerical inputs and
plot specifications are Julia values; a small C++ bridge runs ROOT in isolated
batch processes and writes PNG/PDF. No Julia plotting package, Python, PyCall, or
GUI is required. ROOT itself must be installed separately.

**Status:** experimental, local development only; not published or registered.
The former package name was `JuliaROOT`. See [migration](docs/migration.md).

![Synthetic ROOT graphics gallery](docs/images/overview.png)

*Native ROOT output from `scripts/graphics_gallery.jl`; synthetic inputs only.*

## Requirements and installation

- Julia **1.11+**.
- CERN ROOT **6.40+**, including `root-config`, headers and libraries.
- A ROOT-compatible C++ compiler; `c++` by default, or the executable set in `CXX`.
- macOS/Linux are the development targets. Windows is not validated.

From this source checkout:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/doctor.jl
```

From a separate Julia analysis environment:

```julia
using Pkg
Pkg.develop(path="/absolute/path/to/ROOTPlots-checkout")
using ROOTPlots
```

This is a path dependency, not `Pkg.add("ROOTPlots")` from General. Detailed
setup, cache behavior, and server troubleshooting: [installation](docs/installation.md).

## Quick start

```julia
using ROOTPlots

h = histogram([.2, .4, 1.2, 1.7], [0., 1., 2.]).hist
files = render_root(h, GraphSpec(xlabel="x", ylabel="Counts"),
                    "output/histogram"; draw=:e1, format=:both)
@show files.png files.pdf

observed = GraphSeries([0., 1., 2., 3.], [1.1, 2.9, 5.2, 6.8];
                       yerrors=fill(.2, 4), label="Synthetic data")
fit = fit_root(observed, :pol1; initial=[0., 1.])
fit.status == 0 && fit.valid || error("Fit did not converge")
render_root([fit_curve(fit), observed], GraphSpec(xlabel="x", ylabel="Response"),
            "output/fit"; format=:both)
```

`format=:png` is the default. `:pdf` and `:both` are also supported; unselected
return paths are `nothing`. All example data are **synthetic, not measurements**.

## Scope

| Task | Main API |
| --- | --- |
| Weighted 1D/2D/3D histograms, flows, projection | `histogram`, `Hist1D/2D/3D`, `projection` |
| Points, curves, asymmetric errors and bands | `GraphSeries`, `GraphSpec`, `GraphStyle` |
| Heatmaps, contours, surfaces, XYZ/error graphs | `Grid2D/3D`, `Graph3D`, `Graph3DErrors` |
| Stacks, profiles, convex polygon bins | `HistStack`, `Profile1D/2D/3D`, `PolygonHist2D` |
| ROOT fits, covariance and residual diagnostics | `fit_root`, `FitResult`, `chi2_pvalue` |
| Native TF1, spline, QQ, polar, efficiency, pie | `ModelCurve`, `SplineCurve`, `QQGraph`, etc. |
| Mixed panels, annotations, axes, legends | `PlotFigure`, `PlotPanel`, `AxisStyle`, `PlotText` |
| Data/model comparisons with ratio panels | `comparison`, `Component`, `PanelSpec` |
| Bounded parallel independent ROOT jobs | `rootmap` |

Not a complete ROOT binding. No ROOT-file I/O, RDataFrame, GUI/OpenGL, RooFit,
RooStats, or arbitrary C++/TF1-expression evaluation. Use tools such as UnROOT
for input files. Nested layouts, arbitrary analytic TF2/TF3, nonconvex polygons,
and multi/bent-error graphs are not implemented. See the
[API and scientific semantics](docs/api.md) for supported options and limits.

## Dependencies

**CSV.jl and DataFrames.jl are included** for convenient input and analysis.
Use vectors directly or `histogram(df, :x, edges; weight=:w)` and the 2D/3D column
adapters. No extension setup or separate examples environment is needed.

SpecialFunctions.jl supplies error/incomplete-gamma functions; Scratch.jl
manages a writable native build cache outside the package source. SHA and
Printf/Statistics are standard-library dependencies. No Makie, CairoMakie, or Plots
dependency is present, even optionally.

## Examples

Examples (native ROOT required):

```sh
julia --project=. scripts/demo.jl
julia --project=. scripts/graphfit_demo.jl
julia --threads=4 --project=. scripts/graphics_gallery.jl
julia --project=. scripts/special_graph_demo.jl
julia --project=. scripts/edge_profile_demo.jl
```

Outputs go to `examples/output/`. CSV/DataFrames analysis examples use the same
project environment:

```sh
julia --project=. scripts/histogram_fit_demo.jl
```

## Validation

```sh
# Fails if a Julia plotting package is in the dependency graph; checks CSV/table input.
julia --project=. test/dependency_smoke.jl
# Julia numerical/input tests; native ROOT sections are explicitly skipped.
ROOTPLOTS_TEST_ROOT=off julia --threads=4 --project=. -e 'using Pkg; Pkg.test()'
# Full integration: fails, rather than silently skips, if ROOT is absent.
ROOTPLOTS_TEST_ROOT=required julia --threads=4 --project=. -e 'using Pkg; Pkg.test()'
# Visual regression against reference PNGs (see docs/testing.md).
ROOTPLOTS_VISUAL=required julia --project=. -e 'using Pkg; Pkg.test()'
```

Tests cover numerical results, invalid inputs, native output, concurrent
builds/jobs, and pixel comparison of representative figures against reference
PNGs. Pixel references are only comparable in the environment that produced them,
and passing tests are not proof of scientific equivalence. See
[testing](docs/testing.md) for modes and validated environments.
The included GitHub Actions workflow is a configuration for future use; no
remote CI result is implied before it actually runs.

## License

[MIT](LICENSE) for the original package code. ROOT and other external
dependencies retain their own licenses; see [third-party notes](THIRD_PARTY.md).
