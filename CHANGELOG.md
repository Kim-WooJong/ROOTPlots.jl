# Changelog

## 0.2.1 (2026-10-30)

- Design principle recorded (README, `docs/design.md`, API docs): ROOTPlots is a Julia
  API for ROOT's own drawing and fitting; every option maps onto a ROOT
  method/option/style setting and is named after it; unset options keep ROOT's
  default behaviour.
- Experiment styles, opt-in: `style=:Name` on `GraphSpec`, `Graph2DSpec`, `Graph3DSpec`,
  `PanelSpec` and `PlotFigure` selects a ROOT `TStyle` with `gROOT->SetStyle` before any
  object is created (ROOT's built-ins such as `:ATLAS`, `:Plain`, `:Modern`, and a `"CMS"`
  style built in the backend from ROOT setters following the public CMS
  publication style conventions). Names are checked against ROOT's style list
  (`root_styles()` lists them). Explicit axis/canvas/legend/stats options still win.
- `cms_label`: "CMS" (font 61) with extra text in italics (font 52, default
  "Preliminary") inside the frame at the top left (or `:inside_right`, `:outside`), and
  the luminosity/energy text above the frame on the right ("138 fb^{-1} (13 TeV)").
  `style=:CMS` draws it by default (`cms_label=false` removes it). `atlas_label`:
  "ATLAS" (font 72) followed by "Internal"/"Preliminary" (font 42) and an optional
  "#sqrt{s}, luminosity" line. Labels are ROOT TLatex, sized per pad; they work on
  1-D graphs/histograms, Hist2D COLZ, comparison ratio panels (the label carries
  `lumi`/`annotation`), fit/pull views (upper pad) and every `PlotFigure` pad.
- With a style, two-pad views (fit/pull, comparison) take the style's margins,
  ticks and text sizes rescaled to each pad's height; flat 2-D views keep the
  style's margins (palette room is kept). Without a style nothing changes; the
  existing visual references are unchanged.
- Palette room follows the z text sizes: a COLZ palette's right margin and z title
  offset are computed from the z label/title sizes, the widest palette label and the
  `×10^n` exponent (measured with ROOT's TLatex), after explicit z-axis and canvas
  styles. Under `style=:CMS` (title 0.06) the z title was pushed past the right edge
  of the canvas or figure pad; under ATLAS (top margin 0.05) the exponent was cut at
  the top. ROOT's default style keeps the 0.2 margin and 1.3 offset while the
  labels fit, and explicit `CanvasStyle` margins and z `title_offset` still win.
  `ROOTPLOTS_AXIS_DUMP` records `PALETTE` rows (palette box, z-title extent) and each
  pad's place in the canvas; tests check them and the painted pixels.
- Experiment labels above the frame fit inside their pad. Under `style=:ATLAS` (top
  margin 0.05) the outside ATLAS label (the default on 2-D views) reached the top
  edge of the canvas; the "#sqrt{s}" line was cut. Texts above the frame are now
  measured with ROOT's TLatex and, unless `CanvasStyle(top_margin=...)` is given,
  the pad's top margin grows before anything is painted until they end a small
  white edge below the pad top (single canvases, figure pads, fit/pull and
  comparison upper pads; a later palette or side-margin change re-checks it). An
  explicit top margin is kept and the texts move down instead. The outside ATLAS
  row stacks its second line under "ATLAS <text>" when the pad is too narrow for
  both on one row. The ATLAS baseline sits a quarter text size above the frame.
  CMS labels already fit (their sizes follow the top margin) and are unchanged;
  the visual references are unchanged. `ROOTPLOTS_AXIS_DUMP` `LATEX` rows gain the
  TLatex box width and height; new tests check the dump and the painted pixels for
  ATLAS and CMS labels on every view.
- `ROOTPLOTS_AXIS_DUMP` records `GSTYLE` and `LATEX` rows; new tests
  `test/experiment_style_tests.jl`, a `cms_style` visual reference, and
  `scripts/experiment_style_demo.jl`.

- ROOT axis number/label options in `AxisStyle`, each one ROOT setter:
  `no_exponent` (`TAxis::SetNoExponent`), `more_log_labels` (`SetMoreLogLabels`),
  `decimals` (`SetDecimals`), `optimize_divisions` (`SetNdivisions(n, optim)`),
  `center_labels`, `label_color`/`title_color`/`axis_color`, `ticks` (`SetTicks`),
  `time_display`/`time_format`/`time_offset` (time axes), `exponent_offset`
  (`TGaxis::SetExponentOffset`), `bin_labels` (`SetBinLabel`) with `labels_option`
  (`LabelsOption`), and `change_labels` (`ChangeLabel`, new `LabelChange` type).
  Unset fields keep ROOT's defaults.
- `CanvasStyle`: `grid_x`/`grid_y` (`TPad::SetGridx/y`) and an opt-in ROOT title box
  (`opt_title=true`, with `title_x/y/w/h`, `title_align`, `title_font`,
  `title_font_size`, `title_border_size`, `title_fill_color`, `title_style`,
  `title_text_color` mapping to the gStyle title setters). The default TLatex
  heading is unchanged; ROOT's "title;x;y" form works in ROOT-title mode.
- `Graph2DSpec`/`Graph3DSpec` gain `xaxis`, `yaxis`, `zaxis` (also the palette axis),
  and `canvas`, applied to Grid2D/Grid3D, Hist2D/Hist3D, profiles, TH2Poly, XYZ
  graphs, and TF2/TF3. `PanelSpec` gains `xaxis`, `yaxis`, `ratio_yaxis`, and
  `canvas` (grids, ticks); the fit/pull figure now accepts axis styles, grids,
  ticks, and the title box. All of it works per pad in `PlotFigure`.
- Test hook `ROOTPLOTS_AXIS_DUMP` records ROOT's pad/axis state; new demo
  `scripts/axis_format_demo.jl`; new tests `test/axis_format_tests.jl`.
- Polar and pie charts now reject axis styles, grids, and the title box instead of
  silently ignoring them.
- ROOT statistics boxes: `stats=true` shows ROOT's default box (name, Entries,
  Mean, Std Dev); fields can be chosen as a tuple
  (`stats=(:entries, :mean_error, :underflow, :overflow, …)`), as ROOT letters
  (`"ksiourmen"`), or as a ROOT OptStat integer. Works for every `Hist1D` draw
  mode (presets included), Hist1D overlays (`SAMES` layers get stacked boxes),
  `Hist2D`/`Hist3D`, profiles, `PolygonHist2D`, and panels in `PlotFigure`.
  Previously only raw-option Hist1D/profile boxes without Entries were available;
  `stats=true` now shows ROOT's default content instead of Mean/Std Dev only.
- ROOT fit boxes: `attached_fit=FitResult` attaches a native TF1 with the fitted
  parameters, errors, χ², and ndf; `fitstats=true` (or fields/`"pcev"`/integer)
  shows χ²/ndf, Prob, and parameters ± errors with ROOT's parameter names.
  Available for Hist1D (all draw modes), graphs (`render_root(series, …;
  attached_fit, attached_to, fitstats)`), and figure panels; `fit_label` adds the
  TF1 to the legend.
- Box placement: `stats_position` (`:auto` = ROOT's default, corners inside the
  frame, or an NDC box) and `stats_name`; automatic legends move below the box;
  boxes are widened when ROOT would shrink long lines; 2D boxes avoid the palette.
- `histogram` now records ROOT's Entries (all fills, including flows), the flow
  bins, and the unbinned in-range moments in the returned `Hist1D`/`Hist2D`/
  `Hist3D` (new optional fields `entries`, `underflow`/`overflow` or
  `bins_with_flow`, `moments`); native drawing fills ROOT's under/overflow bins.
  Histograms built from bin contents keep working unchanged (ROOT then shows
  effective entries and bin-centre moments). `rebin` keeps these fields.
- `render_root(layers::AbstractVector, spec, stem; kwargs...)` forwards panel
  keywords; unknown panel keywords are named in the error message.
- Test hook `ROOTPLOTS_STATS_DUMP` records the drawn boxes; new demo
  `scripts/stats_box_demo.jl`.
- `stats_headroom=true` (opt-in; default keeps ROOT's behaviour) raises the
  automatic y maximum of 1-D histograms, profiles, graphs, and 1-D panels so
  data stay below a top statistics/fit box.
- Automatic legends use a fixed ROOT-like text size (0.035 of the pad, 0.045 in
  the upper pad of the fit/pull figure) with the box sized from the entries.
  Previously TLegend scaled its text to fill a fixed box, so short labels were
  drawn about twice the axis label size and long labels much smaller.
- `Profile1D` error bars and markers are drawn in one colour (black). ROOT's
  default profile attributes (blue line, black marker) made batch PNG output
  alternate the error-bar colour bin by bin.
- Titled colour palettes (`COLZ` etc. with a z title) get a 0.2 right margin
  and a z-title offset of 1.3, so the title no longer touches the palette
  labels or the pad edge; explicit `CanvasStyle` margins and z-axis styles
  still win. Figure panels now pick their right margin natively.

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
