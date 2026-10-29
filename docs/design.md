# Design

[Quick start](../README.md) · [API](api.md) · [Installation](installation.md) · [Testing](testing.md)

## Principle

ROOTPlots.jl is a Julia API for using ROOT's own drawing and fitting features,
not a reimplementation. Every option maps onto a ROOT method/option/style
setting, named so a ROOT user recognises it (snake_case of the ROOT method;
docstring cites it). An unset option keeps ROOT's default behaviour. Drawing,
statistics, fitting and formatting are done by ROOT; Julia only validates input
and passes data/options.

## What this means in practice

- **Names.** A keyword is the snake_case form of the ROOT call it makes, and its
  docstring names that call: `no_exponent` is `TAxis::SetNoExponent`,
  `more_log_labels` is `TAxis::SetMoreLogLabels`, `grid_x` is `TPad::SetGridx`,
  `stats_position` places ROOT's `TPaveStats`, `style=:ATLAS` is
  `gROOT->SetStyle("ATLAS")`.
- **Defaults.** Option fields default to `nothing`, which means "do not call the
  setter": ROOT's default (or the selected ROOT `TStyle`) stays in place. Where
  the package needs a layout of its own (two-pad ratio and pull views, palette
  room), the documentation says so, and explicit options still win.
- **ROOT does the work.** Axes, labels, number formats, statistics boxes, fit
  boxes, fits, palettes, legends and text are produced by ROOT objects in a
  native batch process (`TH1`, `TGraph`, `TF1`, `TPaveStats`, `TLegend`, `TLatex`,
  `TStyle`, ...). Julia does not draw, and does not recompute what ROOT reports.
- **Julia validates and passes data.** Julia checks inputs early, with clear
  errors, and sends numbers and options to the native bridge in a small text
  protocol. Numeric inputs (histogram contents, the bins and ratios of a
  data/model comparison) are Julia data that ROOT then draws.
- **Opt-in.** New features are off unless requested, so existing plots and the
  visual references do not change.

## Adding an option

1. Find the ROOT method, draw option, or style setter that provides the feature.
2. Name the keyword after it (snake_case) and cite it in the docstring.
3. Default to `nothing` (not set); validate in Julia; pass the value through the
   request; apply it with the ROOT call in the backend.
4. Test that the value reaches ROOT (the `ROOTPLOTS_AXIS_DUMP` and
   `ROOTPLOTS_STATS_DUMP` hooks record ROOT's own state) and that unset options
   leave the output unchanged.
