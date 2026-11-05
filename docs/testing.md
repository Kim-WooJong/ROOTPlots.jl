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

## Native request validation

Malformed requests sent directly to a native backend (both `rootplots-graphfit`
and the comparison renderer `rootplots-render`) end with a controlled error:
exit code 2, no signal, and no PNG or PDF written. The test helpers
`run_native` and `rejects_request` (`test/setup.jl`) run the executable with a
bounded timeout and check exactly that. `test/native_validation_tests.jl` covers
the comparison protocol (missing panels, row order and indices, duplicates,
field counts, non-finite and out-of-range values, negative components, unknown
rows) and the `CMSCANVAS` rows (only directly after a scene's CMS `ROOTSTYLE`
row; field count and value ranges).

The optional `scripts/fuzz_native_requests.jl [cases] [dir]` (not part of
`Pkg.test`) mutates a valid comparison request and a CMS graph request
(deleted, duplicated or swapped rows; fields replaced by empty, non-numeric,
non-finite, negative or huge values) and requires every run to either succeed
with an image or exit with code 2 without one; failing requests are kept.

## Text checks without pixels

A changed axis title, legend label, or luminosity text can stay within the
pixel tolerance of the visual references. `test/text_semantic_tests.jl` therefore
renders the visual-reference figures with the `ROOTPLOTS_AXIS_DUMP` hook and
compares the texts ROOT holds (axis titles, TLatex texts, legend entries and
options) with explicit expectations, and checks statistics-box lines through
`ROOTPLOTS_STATS_DUMP`. No fonts or geometry are compared, so these checks run
wherever native ROOT runs, including platforms whose pixel comparison is skipped.

## Statistical cross-checks

`test/histogram_ops_tests.jl` compares scaling, merging and flow-aware
projections with direct fills of the equivalent events (signed and zero weights,
flow bins) and with numbers printed by plain ROOT 6.40.04 macros.
`test/fit_error_tests.jl` compares ROOT's fitted parameters, covariance and χ²
with closed-form weighted least squares (pol1: uniform and non-uniform errors,
no errors, histogram bin centres with zero-error bins, ignored x errors) and
with the analytic Jacobian `inv(J'WJ)` of a noiseless Gaussian (relative
tolerance 1e-3 for Minuit's finite-difference Hessian).

## Native processes, formats and records

`test/runtime_tests.jl`: without ROOT (POSIX shell fixtures), timeouts within bounded
time (limits of 1 to 3 s, generous upper bounds), the child's own process group,
grandchild termination, SIGKILL after an ignored SIGTERM, 3 MB of stderr while
timing out, disabled limits, interrupts (also while the output streams of an
exited child are drained: no descendant survives), no leaked descriptors, stray
stream holders, removal of a staged record file when its writer fails, per-batch `rootmap` limits (one hung job among successful ones gives an
indexed `RootMapError`), cancellation isolation between concurrent and nested
batches, and the `root-config` metadata cache with a fake `root-config` call counter
(hit, alias path, edited file, force, concurrent first use) and backend keys
(source/header edits, standard, compiler). With ROOT: `--formats` on both
executables (only the requested files, legacy form writes both, PNG
byte-identical to the two-format PNG, PDF equal apart from ROOT's dates/title,
bad options exit 2 without output), a PDF `/Title` of `rendered.pdf` without any
directory for graph, comparison, figure and histogram routes, and records for graph, histogram, step
histogram, comparison, nested figure, fit/pull and TH2 routes (TOML parses,
digests match, replay gives a byte-identical PNG, no directory/user/host marker,
no record by default or on failure).

## Validated environments

| Date | OS | Julia | ROOT | Result |
| --- | --- | --- | --- | --- |
| 2026-11-05 | macOS (arm64) | 1.12.6 | 6.40.04 | full native suite with visual comparison required, after the histogram operations and validation work (2028/2028; references unchanged); rootless 1051 pass + 10 explicit skips; request fuzzing 1607 mutated requests without failure |
| 2026-11-05 | macOS (arm64) | 1.12.6 | 6.40.04 | full native suite with visual comparison required, after the correctness fixes (1665/1665; references unchanged); rootless 879 pass + 3 explicit skips |
| 2026-11-04 | macOS (arm64) | 1.12.6 | 6.40.04 | full native suite with visual comparison required, after the cmsstyle alignment (1507/1507; references cms_style regenerated, cms_outside added); rootless 780 pass + 2 explicit skips |
| 2026-10-29 | macOS (arm64) | 1.12.6 | 6.40.04 | full native suite incl. experiment styles (1350/1350), visual references recorded here |
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
