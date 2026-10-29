# Migrating JuliaROOT to ROOTPlots (0.2)

This is a deliberate package rename and dependency reduction, made before the
first public release, not a compatibility alias. The package UUID is unchanged.

1. In an analysis environment that tracks the old name, remove that dependency
   before adding the renamed package (from GitHub, or a local checkout):

   ```julia
   using Pkg
   Pkg.rm("JuliaROOT")  # only where the old name is already a dependency
   Pkg.add(url="https://github.com/Kim-WooJong/ROOTPlots.jl")
   # or: Pkg.develop(path="/absolute/path/to/checkout")
   Pkg.instantiate()
   ```

2. Change `using JuliaROOT` and qualified `JuliaROOT.` references to `ROOTPlots`.
   Restart Julia after changing the environment/name.
3. `render_cairo` and Cairo fallback are removed. Use `render_root`; all rendering
   requires native ROOT. Numeric input/binning operations still work without it.
4. CSV and DataFrames remain included for input and analysis. DataFrame and vector
   histogram APIs remain supported. All repository scripts use `--project=.`;
   no optional extension activation or separate examples environment is needed.
5. `scripts/demo.jl [output-stem]` is now ROOT-only. The old `auto|root|cairo`
   backend argument is removed.
6. Native artifacts now use a depot scratch cache. Old checkout-local `build/`
   files are no longer used; they are not automatically deleted. Rebuild after
   installing/updating the sources. Private TSV protocol names changed together
   with both bridge executables; do not mix old/new backend files.
7. Behaviour changes within 0.2.0: `PolygonHist2D` now accepts simple
   non-convex polygons and drops a repeated closing vertex (both were rejected
   before); self-intersecting polygons are still rejected. `PlotFigure.panels`
   may contain nested `PlotFigure`s (element type `Union{PlotPanel,PlotFigure}`).
   ROOT `Info` messages are no longer printed; set `ROOTPLOTS_ROOT_VERBOSE=1`
   to see them.

The development checkout directory is `ROOTPlots.jl`. Existing copies elsewhere
(for example an older checkout named `julia-ROOT`) keep their old path and version
until explicitly synchronized and migrated. Package identity comes from
`Project.toml` and `src/ROOTPlots.jl`, not the directory spelling.
