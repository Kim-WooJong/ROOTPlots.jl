# Installation and runtime

## Choose an existing ROOT installation

ROOTPlots builds two small C++ executables against your ROOT installation; it
does not build or download ROOT. A compatible compiler and ROOT headers/libraries
are required. No X11 connection is needed for batch PNG/PDF output.

```sh
julia --version
root-config --version
root-config --cflags --libs
c++ --version
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/doctor.jl
```

`doctor.jl` checks discovery, compilation/linking, an actual native fit, and PNG/PDF
creation. Read compiler/loader errors from the stage that failed. It is an
environment smoke check, not a validation of your scientific model.

Resolution order is an explicit `root_config` keyword, `ROOT_CONFIG`,
`ROOTSYS/bin/root-config`, then PATH. An invalid explicit selection or
`ROOT_CONFIG` does **not** silently fall back to a different installation.

```sh
ROOT_CONFIG=/absolute/ROOT/bin/root-config julia --project=. scripts/doctor.jl
CXX=/absolute/compiler/bin/clang++ julia --project=. scripts/doctor.jl
```

`CXX` must be a compiler executable, not a shell command or flags. ROOT's own
compile/link flags are obtained from `root-config --cflags --libs`.

```julia
using ROOTPlots
built = build_root()  # also built automatically on first native use
build_root(root_config="/absolute/ROOT/bin/root-config", force=true)
```

## Cache and read-only package sources

Native artifacts live in Scratch.jl's `native-backends` space under the writable
Julia depot (normally `~/.julia/scratchspaces/`), **not** in the package checkout.
`build_root()` returns the actual paths. Source files may be read-only.

The key includes source/header content, ROOT path/version, compiler path,
language standard and ROOT flags. Builds use a process-local lock and unique
staging files, followed by atomic publication. Independent Julia processes may
compile redundantly but do not expose partially written executables. Failed
forced builds leave the previous executable intact.

Use `force=true` after replacing an SDK/compiler or ROOT in place without
changing its path/version. Individual old build keys are not automatically
evicted by ROOTPlots. Scratch/Pkg manage the cache lifecycle; do not delete
artifacts while jobs are running. No source cleanup is needed after building.

## macOS

Install ROOT using the [ROOT installation instructions](https://root.cern/install/)
(Homebrew is one option) and install compatible Command Line Tools. For a
Homebrew installation, check `brew --prefix root` if `root-config` is not on PATH.
ROOT/Cling can depend on SDK paths from its build environment: SDK errors should
be resolved by aligning ROOT, macOS and the developer tools, not by modifying
the package's scientific code.

## Linux / LCG environments

Use ROOT and a compiler from the same compatible stack. A bridge can compile
successfully yet fail at runtime if another `libstdc++`, TBB, VDT, or ROOT library
is selected. Run `ldd` on the paths returned by `build_root()` when diagnosing
missing libraries or `GLIBCXX_*` errors. Avoid mixing system and LCG toolchains.

Setting a stack's `LD_LIBRARY_PATH` before starting Julia also affects Julia and
its own libraries. If the stack requires conflicting libraries, start Julia
with a clean loader environment and set the required library search path in
Julia's `ENV` before launching the ROOT child. Do not change `ENV` while parallel
jobs are active. This package does not install or alter system libraries.

Existing server wrappers can remain in place; updating the sources does not
automatically update other copies or analysis environments. See [migration](migration.md).

## CSV and DataFrames input

CSV and DataFrames are included in the package environment. No separate
extension or examples project needs to be set up:

```julia
using ROOTPlots, DataFrames, CSV
h = histogram(DataFrame(x=[.2, 1.2]), :x, [0., 1., 2.]).hist
```

For repository examples that read CSV, use `--project=.` like the other scripts.
In a separate analysis project, explicitly add CSV/DataFrames to that project's
dependencies if your own code imports them directly; ROOTPlots' transitive
dependencies are not a substitute for declaring your analysis dependencies.

The current support target is Julia 1.11+ / ROOT 6.40+ on macOS/Linux. A version
range is not proof that every combination works. See [testing](testing.md)
for the environments actually exercised; Windows, older ROOT, and ROOT 7 are not validated.
