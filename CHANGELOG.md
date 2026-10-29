# Changelog

## 0.2.0 — unreleased, local development

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

This version has not been published. See `docs/migration.md` for breaking changes.
