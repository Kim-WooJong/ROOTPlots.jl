# Third-party dependencies

The MIT license in `LICENSE` covers the original ROOTPlots source and examples,
including its Julia API and C++ bridge. It does not relicense external libraries.
No ROOT distribution or precompiled ROOT libraries are bundled here.

| Dependency | Use | Upstream license/source |
| --- | --- | --- |
| CERN ROOT | External, user-installed plotting and fitting runtime | [LGPL-2.1 and component notices](https://github.com/root-project/root/blob/master/LICENSE) |
| SpecialFunctions.jl | Error functions and incomplete gamma for model evaluation and statistics | [MIT and third-party notices](https://github.com/JuliaMath/SpecialFunctions.jl) |
| Scratch.jl | Writable native build cache outside the package source | [MIT](https://github.com/JuliaPackaging/Scratch.jl) |
| DataFrames.jl | Table-to-histogram adapters and analysis | [MIT](https://github.com/JuliaData/DataFrames.jl) |
| CSV.jl | Read/write input tables | [MIT](https://github.com/JuliaData/CSV.jl) |

Julia's standard libraries and transitive dependencies retain their own licenses.
When redistributing a ROOT runtime or compiled binaries, check the corresponding
upstream component and redistribution requirements separately. Source-only
ROOTPlots does not install or redistribute ROOT.
