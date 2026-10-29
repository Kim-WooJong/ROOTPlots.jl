# Changelog

## 0.2.0

- Rename package/module `JuliaROOT` to `ROOTPlots`; preserve UUID.
- Remove Cairo rendering and all Julia plotting-package dependencies.
- Retain CSV and DataFrames for input/analysis; share vector-based 1D/2D/3D
  numerical kernels in the DataFrame adapters.
- Move native build artifacts outside the package source using Scratch.jl.
- Support a compiler executable in `CXX`; reject invalid explicit ROOT selections.
- Add native-required/rootless test modes, dependency smoke checks, environment
  diagnostics, and CI configuration (not yet remotely executed).
- Add MIT licensing, third-party notes, and consistent installation/API/migration docs.
- English-only documentation; generic synthetic examples (`edge_profile_demo.jl`)
  replace lab-specific scripts; comparison panels no longer add an experiment label.
- Add visual regression tests comparing representative figures with reference
  PNGs through ROOT's image reader (`ROOTPLOTS_VISUAL`).
- Install from GitHub with `Pkg.add(url="https://github.com/Kim-WooJong/ROOTPlots.jl")`
  (not registered in General); a source checkout with `Pkg.develop` remains
  available for development.
- Quiet native ROOT output: ROOT's message level is set to warning, so
  `Info in <TCanvas::Print>` lines no longer appear for every saved file;
  warnings and errors are still printed. `ROOTPLOTS_ROOT_VERBOSE=1` restores
  ROOT's default verbosity.
- Add `BentErrorSeries` (native TGraphBentErrors: asymmetric errors with
  shifted error-bar end points).
- Add `MultiErrorSeries`/`YErrorSet` (native TGraphMultiErrors: several y-error
  sets with per-set display and style, ROOT sum modes, legend entries per set)
  and `total_yerrors`.
- `render_root(layers::AbstractVector, ::GraphSpec, stem)` draws mixed XY layers
  (graphs, bent/multi-error graphs, TF1 curves, histograms) on one pad.
- Add `ModelSurface` (native TF2) and `ModelVolume` (native TF3) from predefined
  named formulas, following the `ModelCurve` design; TF2 supports COL/CONT/SURF/LEGO
  options, TF3 the batch-renderable iso-surface and sampled BOX2Z views. Add
  `function_grid` to sample Julia callables of two or three variables on a grid.
- `PolygonHist2D` accepts simple non-convex polygons (validated for finiteness,
  area, and self-intersection; a repeated closing vertex is dropped). Add the
  `PlotPolygon` filled-polygon primitive.
- `PlotFigure` cells may be nested `PlotFigure`s (genuinely nested TPads, up to
  4 levels) and accept relative `widths`/`heights`; regular figures produce the
  same native request as before.

Published on GitHub; not registered in the General registry. See
`docs/migration.md` for breaking changes.
