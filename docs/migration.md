# Migrating JuliaROOT to ROOTPlots (0.2 development)

This is a deliberate pre-publication package rename and dependency reduction,
not a compatibility alias. The package UUID is unchanged. No remote repository,
registry entry, or release was created by this change.

1. In an analysis environment that tracks the old name, remove that dependency
   before developing the renamed package at its local path:

   ```julia
   using Pkg
   Pkg.rm("JuliaROOT")  # only where the old name is already a dependency
   Pkg.develop(path="/absolute/path/to/checkout")
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

The development checkout directory is `ROOTPlots.jl`. Existing copies elsewhere
(for example an older checkout named `julia-ROOT`) keep their old path and version
until explicitly synchronized and migrated. Package identity comes from
`Project.toml` and `src/ROOTPlots.jl`, not the directory spelling.
