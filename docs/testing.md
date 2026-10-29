# Testing

## Test modes

| Command | Purpose |
| --- | --- |
| `julia --project=. test/dependency_smoke.jl` | Fails if a Julia plotting package is in the dependency graph; checks CSV/DataFrame input. |
| `ROOTPLOTS_TEST_ROOT=off ... Pkg.test()` | Numerical and input-validation tests only; native sections are explicitly skipped. |
| `ROOTPLOTS_TEST_ROOT=required ... Pkg.test()` | Full integration; fails instead of skipping when ROOT is absent. |
| `julia --project=. test/readonly_install.jl` | Fresh analysis environment with a read-only package source. |
| `julia --project=. scripts/doctor.jl` | Environment check: compile, link, native fit, PNG and PDF output. |

`ROOTPLOTS_TEST_ROOT` accepts `auto` (default), `required`, or `off`.

## Visual regression

`test/visual_tests.jl` renders representative figures (TH1 errors, graph band
with TF1, TH2 heatmap, THStack, data/model ratio panel, mixed 2×2 figure) and
compares them with `test/references/*.png` using ROOT's own image reader. A
figure passes when the image size is identical, at most 0.2% of pixels differ by
more than 32/255 in any RGB channel, and the mean channel difference is at most
0.5. A built-in self-check verifies that the comparator detects a real change.

`ROOTPLOTS_VISUAL` selects the mode:

- `auto` (default): compare only when the OS and ROOT version match
  `test/references/environment.txt`; otherwise the comparison is skipped.
- `required`: always compare and fail on missing or different references.
- `update`: re-render and overwrite the references. Review the images before
  committing them.
- `off`: do not run.

Fonts, anti-aliasing, and ROOT painter details differ between platforms and
releases, so references are environment-specific by design.

## Validated environments

| Date | OS | Julia | ROOT | Result |
| --- | --- | --- | --- | --- |
| 2026-10-29 | macOS (arm64) | 1.12.6 | 6.40.04 | full native suite incl. graph extensions (978/978), visual references recorded here |
| 2026-10-29 | Linux EL9 (x86_64), LCG_110 clang 22 | 1.12.2 | 6.40.02 | doctor PASS; full native suite 818/818, visual comparison skipped (different reference environment) |

| 2026-10-29 | GitHub Actions, ubuntu-latest | 1.11, 1 (rootless); 1 (native) | conda-forge 6.40 | all 3 CI jobs passed (before the bent/multi-error, TF2/TF3, polygon, and nested-layout additions) |

Other combinations within Julia 1.11+ / ROOT 6.40+ are supported targets but are
not claimed as validated until they appear in this table. The workflow is
`.github/workflows/ci.yml`.

The graph-extension tests (`test/graph_feature_tests.jl`) cover input
validation without ROOT and, with ROOT, native rendering of every new display
mode, a pixel check that a non-convex polygon is not filled like its convex hull,
malformed-request rejection by the native bridge, and the default absence of
`Info in <TCanvas::Print>` messages (and their return with
`ROOTPLOTS_ROOT_VERBOSE=1`).
