# Changelog

## 0.3.0 (2026-11-05)

- Native process limits: every native process (ROOT backends, `root-config`, the
  compiler) runs under a time limit (`ROOTPLOTS_BACKEND_TIMEOUT` 600 s,
  `ROOTPLOTS_ROOT_CONFIG_TIMEOUT` 60 s, `ROOTPLOTS_COMPILE_TIMEOUT` 600 s; `0`/`none`
  disables; invalid values fail before a process starts). On macOS/Linux the
  process leads its own process group; on a timeout or an interrupt of the waiting
  task the group gets SIGTERM, then SIGKILL after 2 s, and the child is reaped
  (Windows: direct child only). `RootBackendError` gained `timedout` and `timeout`
  fields (the six-field constructor still works). `rootmap(...; timeout=t)` sets a
  per-batch backend limit through a scoped value; an interrupt cancels only that
  invocation's processes (nested batches included), stops scheduling, joins the
  workers and rethrows the `InterruptException`. Julia callbacks are not timed out.
  An interrupt or error at any point of a native run, including the wait for
  output streams still held by descendants after the child exited, stops the
  owned process group and waits for it before the exception is rethrown.
- Compatibility: compiler and `root-config` failures are now `RootBackendError`
  (was `ProcessFailedException`); native processes get `devnull` as standard input
  and run in their own session (a terminal Ctrl-C reaches Julia, which stops them).
- Output formats: `format=:png` / `:pdf` make ROOT write only that file (the native
  renderers take `--formats=png|pdf|png,pdf`; without it both files are written as
  before). Images are unchanged: the PNG is byte-identical to the PNG of a `:both`
  run; visual references unchanged. The renderer now runs in its private temporary
  directory with a relative output stem, so a PDF's `/Title` is `rendered.pdf`
  instead of an absolute temporary path (pixels unchanged). A relative `TMPDIR`
  works: the private directory, request and working directory are made absolute,
  and the child alone gets `TMPDIR` as an absolute path (Julia's `ENV` unchanged).
- Reproducibility records (opt-in): `provenance=true` on every `render_root` method
  writes `<stem>.request.tsv` (exact request, data verbatim) and
  `<stem>.provenance.toml` (versions, OS/arch, backend key, style, formats,
  SHA-256 digests; no paths, flags, host or user names) after a successful
  rendering and returns `(; png, pdf, request, provenance)`. A failure while
  writing either record removes its temporary file; the two renames are not
  atomic as a pair (documented). Default return values
  and files are unchanged. New standard-library dependency: TOML.
- `root-config` output (version, flags) is cached per `root-config` file identity (path,
  real path, inode, size, mtime) for the session; backend sources are still hashed
  and `CXX` resolved on every lookup; `build_root(force=true)` refreshes. Cache hits no
  longer take the build lock.
- New `scripts/benchmark.jl` (warmed medians/IQR; no speed-up is claimed beyond its
  measurements) and `scripts/runtime_demo.jl`.
- Histogram operations (new): `c * h`, `h * c`, `h / c` scale a `Hist1D`/`Hist2D`/`Hist3D`
  like ROOT's `TH1::Scale` (every bin including flows; `Σw²` by c²; moments scaled,
  entries kept; finite factors only, `h / 0` and Float64 overflow are errors;
  `0 * h` keeps the fill count as `Scale(0)` does). `merge(a, b, more...)` combines
  independent samples with identical edges like `TH1::Add(h, 1)` (bins, flows,
  variances, known entries and moments add; unknown stays unknown; the same
  histogram twice is rejected with advice to scale instead). `projection(...;
  flow=false)` gained `flow=true` for 2D/3D → 1D and 3D → 2D: integrated axes
  include their flow bins, projected axes keep theirs, and the original entries
  and moments are kept only when exact (conservative rule; differences from
  ROOT 6.40 `ProjectionX` documented). A `Hist2D`/`Hist3D` without
  `bins_with_flow` counts as having empty flow bins (documented assumption). The
  default projection is unchanged. Scaling does not check Float64 underflow (tiny
  factors can round variances to zero; documented). All results own new arrays. New `docs/recipes.md` and `scripts/histogram_ops_demo.jl`.
- Native request validation: the comparison renderer (`rootplots-render`) now
  exits with code 2 for malformed requests, like `rootplots-graphfit` (was 1), and
  rejects more malformed input: duplicate `PCANVAS`/`PAXIS` rows, `BIN` rows after
  `COMP` rows or beyond the bin count, `COMP` rows before all bins, invalid
  component labels/colours, negative errors or predictions, and validity flags
  other than 0/1. `rootplots-graphfit` rejects a `CMSCANVAS` row that does not
  directly follow a scene's CMS `ROOTSTYLE` row (it was silently ignored), and
  checks its field count and ranges in every request. Requests written by the
  package are unchanged.
- `AxisStyle(time_offset=t)` requires `abs(t) < 2^31` seconds (Julia and the native
  bridges): ROOT 6.40's `TAxis::SetTimeOffset` keeps 32-bit seconds, printed a wrong
  offset beyond 2038 and crashed (SIGSEGV) on huge values. Found by the new optional
  `scripts/fuzz_native_requests.jl`, which mutates valid requests and requires every
  run to succeed with an image or fail with exit code 2.
- Test hook: `ROOTPLOTS_AXIS_DUMP` also writes `LEGEND` (box) and `LEGENTRY`
  (label, option) records. New tests check the texts of the visual-reference
  figures without pixels, the native exit-code contract, and fit covariances
  against closed-form least squares and the analytic Jacobian.
- The CMS style's upstream reference is pinned to cms-cat/cmsstyle commit
  `ecfed5a18a7a8da5e49bfce903e40e0f5769c4f1` (style values unchanged).
- `scripts/correctness_demo.jl`: the CMS example with the long y title
  "Events / GeV" uses cmsstyle's recipe `cms_label=(extra_space=0.05,
  y_title_offset=1.1)` (with the default non-square geometry the title reached
  the canvas edge); the CMS default is unchanged.
- Log-y histogram fit/pull views no longer fail on ordinary data with empty tail
  bins or single-entry bins. The upper pad now draws the histogram as a native
  ROOT `TH1D` (`E1 X0`: markers and capped error bars, no bin-width x bars)
  instead of a graph, so ROOT's TH1 rules apply: empty bins are not drawn (also on
  linear axes, where a marker at zero was drawn before), bins below a log frame are
  skipped, and error bars reaching zero are clipped at the frame. The automatic log
  range follows the drawn bins; the fit curve is clipped below them. A log axis
  without any positive bin, or with non-positive `ylimits`, is still an
  `ArgumentError`. Protocol: `DRAWFIT` requests carry a `FITHIST` block.
- Graphs on log axes: error bars that reach zero or below are accepted and
  clipped at the frame edge by ROOT (before, `render_root` rejected them); point
  centres must still be positive. The automatic range follows the bars' positive
  part. Applies to `GraphSeries` in graphs, panels and overlays (bent and
  multi-error graphs keep their stricter checks).
- Unicode labels: every text sent to ROOT (titles, axis titles, legend entries,
  statistics-box headers, `PlotText`, line labels, CMS/ATLAS texts, bin labels,
  `LabelChange` texts, comparison-panel texts, time formats) is translated to ROOT
  TLatex: all Greek letters, superscript/subscript digits and signs (`GeV/c²` →
  `GeV/c^{2}`), `± × ÷ · ° − ≤ ≥ ≠ ≈ ∝ ∞ ∂ ∇ √ ∑ ∫ ∈ ⊂ ∩ ∪`, arrows, `ℏ`, `Å`,
  accented letters (`x̄` → `#bar{x}`, `é` → `#acute{e}`) and typographic
  punctuation. Before, only `μ` and `ℓ` were translated (`ℓ` to `#ell`, which
  TLatex does not know; it is now `#font[12]{l}`) and any other Unicode label was
  drawn as mojibake. ASCII text (including ROOT `#` syntax) is unchanged, also
  when mixed with translated characters: `"#sigma_{m} [GeV/c²]"` →
  `"#sigma_{m} [GeV/c^{2}]"`, `"#splitline{α}{β}"` → `"#splitline{#alpha}{#beta}"`
  (only a translated `#in`/`#subset`/`#supset` followed by a letter gets ROOT's
  `^{}` separator). `GraphSeries` stores any label without tabs, newlines or NUL
  (as in 0.2.2); drawing translates it and rejects text TLatex cannot draw, while
  `fit_root` ignores the label (it draws nothing), so series labelled with e.g.
  Hangul or emoji can still be fitted.
  Characters TLatex cannot draw (e.g. Hangul, CJK, emoji) now raise an
  `ArgumentError` naming the character before ROOT runs, and a request file that
  still contains non-ASCII text is refused.
- Native failures: a failed ROOT process now throws the new exported
  `RootBackendError` (operation, request kind, executable, exit code or signal,
  and the process's stderr text, bounded to its first and last 32 KiB) instead of
  `Base.ProcessFailedException` with the message only printed to the terminal.
  Applies to rendering, fitting, profile statistics, special-graph points,
  `root_styles` and the image helpers; `rootmap` keeps each job's
  `RootBackendError` and indents its message under the job index. Standard error
  is read concurrently (large output cannot block a process) and is written to
  Julia's `stderr` after a successful run, so ROOT warnings and
  `ROOTPLOTS_ROOT_VERBOSE=1` output stay visible. Compiler and `root-config` failures also use `RootBackendError`
  (see the compatibility note above).
- `fit_root` rejects more than 1,000,000 histogram bins or graph points in Julia,
  before starting ROOT (the native fit's limit).
- ROOT Entries of a `Hist1D`/`Hist2D`/`Hist3D` with unbinned `moments` but unknown
  `entries`: ROOT showed the number of bins the bridge had set (visible plus flow
  bins). It now shows ROOT's effective entries from `TH1::ResetStats` (as without
  moments; ROOT 6.40 includes the flow bins), while Mean/Std Dev keep the supplied
  unbinned moments. Known entries are unchanged.
- Test hook: `ROOTPLOTS_AXIS_DUMP` also writes `DRAWN` records (class, name, draw
  option) for every pad primitive. New script `scripts/correctness_demo.jl`.

- CMS style now follows the CMS Collaboration's cmsstyle package
  (https://github.com/cms-cat/cmsstyle), the reference for CMS-style figures:
  - `style=:CMS` TStyle values are cmsstyle's `setCMSStyle`: pad margins top 0.05,
    bottom 0.13, left 0.16, right 0.02 (were 0.08/0.13/0.16/0.04), axis title offsets
    x 1.1, y 1.35 (were 0.9/1.25), label offset 0.012 (was 0.007), function line
    width 1 (was 2), marker size 1, legend text size 0.04, hatch width 2 / spacing 1.3.
  - Single-pad canvases (graphs, 1-D/2-D histograms, functions, special graphs,
    one-panel figures) take cmsstyle's `cmsCanvas` geometry: top 0.07, bottom
    0.145, right 0.05 H/W, left from 0.145 H and the y title offset, x title offset
    1.05, y title offset 0.78 (non-square) or 1.2 (square canvas). Two-pad views,
    multi-panel figures and 3-D views keep the `setCMSStyle` margins.
  - The CMS label uses cmsstyle's `CMS_lumi` geometry: inside positions 0.035 of
    the frame from its edges (was 0.045), additional lines at the extra-text size
    and cmsstyle's spacing, the outside extra text at cmsstyle's offset after "CMS".
  - New `cms_label` keys: `right` (free text at the luminosity position),
    `additional_info` (AppendAdditionalInfo), `cms_text` (SetCmsText), `scale_lumi`,
    `extra_space` and `y_title_offset` (cmsCanvas extraSpace/yTitOffset; need
    `style=:CMS`); `position` also accepts `:inside_center` and cmsstyle iPos codes
    0/11/22/33; `extra` accepts the nicknames "p", "s", "su", "wip", "pw", and an
    extra text containing "Private" drops the "CMS" text (cmsstyle behaviour).
  - Protocol: `EXPLABEL` row gained the CMS text and luminosity scale fields; new
    `CMSCANVAS` scene row.
  - Visual reference `cms_style.png` regenerated; new reference `cms_outside.png`
    (ME0-like layout: outside label, right text, extra space, log y). Other references
    unchanged.

## 0.2.2 (2026-10-30)

- Two-pad views (comparison ratio panels and fit/pull), with the default style,
  CMS, ATLAS or any ROOT style: labels at the frame edges between the pads are no
  longer cut. The ratio/pull pad's top margin and the upper pad's bottom margin now
  leave room for half a y label, so e.g. "1.5" and "0.5" of a 0.5-1.5 ratio axis
  and "0" of a linear upper axis are drawn whole (before, the ratio frame touched
  the top of its pad and "1.5" was clipped).
- Two-pad views: one text size in both pads. ROOT sizes axis text relative to each
  pad's height, so the short ratio pad enlarged its labels and titles and the
  upper pad shrank them (most visible when `yaxis` sizes were set: a small
  "Events" title over a large "Data / pred."). The lower pad's labels and titles
  now take the upper y axis's sizes in pixels (under a style, the style's one-pad
  sizes), the ratio/pull y title is placed at the upper title's distance from the
  axis so the two titles line up, and a ratio y title too long for its pad is
  reduced to fit (the former fixed 80% reduction under a style is gone). ROOT's
  automatic y title offset (0 in ROOT's default style, used with a label but no
  style) would place the two titles differently; the upper title then sits next to
  the left pad edge. Explicit `AxisStyle` sizes and offsets still win and keep
  ROOT's meaning for their own pad (documented under "Two-pad views" in the API
  docs; precision-3 fonts give pixel sizes).
- `scripts/experiment_style_demo.jl`: a CMS data/MC comparison with log y and a
  0.5-1.5 ratio panel (`stack_ratio_cms_logy`).
- Visual reference `ratio_panel.png` regenerated (ratio pad text and margins).

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
