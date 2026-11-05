# Time limits, cancellation-safe batches, requested formats, and reproducibility records.
#
#   julia --threads=4 --project=. scripts/runtime_demo.jl [output/directory]
#
# Synthetic data. Part 1 uses a POSIX "sleep" as a stand-in for a native ROOT
# process that hangs (run through the package's internal process runner, the
# same path every ROOT backend takes).
using ROOTPlots, TOML
using ROOTPlots: run_backend, sha256

output = isempty(ARGS) ? joinpath(@__DIR__, "..", "examples", "output", "runtime") : ARGS[1]
mkpath(output)
built = build_root()
digest(path) = bytes2hex(open(sha256, path))

# 1. A batch with one hung native job: the limit stops it, the others finish.
x = collect(0.5:1.0:19.5)
samples = [Hist1D(collect(0.0:1.0:20.0), 50 .* exp.(-((x .- 10) ./ (2 + k)) .^ 2 ./ 2) .+ 1,
                  50 .* exp.(-((x .- 10) ./ (2 + k)) .^ 2 ./ 2) .+ 1) for k in 1:4]
elapsed = @elapsed failure = try
    rootmap(1:4; workers=4, timeout=3) do k
        if k == 2 && Sys.isunix()
            run_backend(Cmd(["/bin/sh", "-c", "echo 'stand-in backend started' >&2; sleep 120"]);
                        operation="stand-in for a hung ROOT job")
        end
        render_root(samples[k], GraphSpec(xlabel="m [GeV]", ylabel="Events"),
                    joinpath(output, "batch-$k"); draw="E1")
    end
    nothing
catch e
    e
end
println("1. rootmap(...; timeout=3) returned after ", round(elapsed; digits=1), " s")
if failure isa RootMapError
    showerror(stdout, failure); println()
    for (i, f) in failure.failures
        println("   job ", i, ": ", typeof(f.ex), ", timedout = ", f.ex.timedout)
    end
end
println("   PNGs of the other jobs: ", filter(f -> startswith(f, "batch-"), readdir(output)))

# 2. Only the requested formats are written.
for format in (:png, :pdf, :both)
    files = render_root(samples[1], GraphSpec(xlabel="m [GeV]", ylabel="Events"),
                        joinpath(output, "format-$format"); draw="E1", format)
    println("2. format=", rpad(repr(format), 6), " -> ", files, "  on disk: ",
            filter(f -> startswith(f, "format-$format."), readdir(output)))
end

# 3. A reproducibility record and its replay.
h = samples[3]
fit = fit_root(h, :gaus; initial=[50.0, 10.0, 4.0])
files = render_root(h, fit, GraphSpec(xlabel="m [GeV]", ylabel="Events", style=:CMS),
                    joinpath(output, "fit-record"); provenance=true)
println("3. ", files)
print(read(files.provenance, String))
record = TOML.parsefile(files.provenance)
exe = record["backend"]["executable"] == "rootplots-render" ? built.renderer : built.graphfit
replay = joinpath(output, "fit-record-replayed")
run(Cmd([exe, files.request, replay, "--formats=png,pdf"]))
println("   replayed PNG identical: ", digest(replay * ".png") == record["outputs"]["png"]["sha256"],
        "; replay also wrote a PDF: ", isfile(replay * ".pdf"))
