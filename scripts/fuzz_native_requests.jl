# Optional request fuzzing of the native backends (not part of Pkg.test).
# Mutates valid requests deterministically (delete, duplicate or swap rows; replace
# one field by an empty, non-numeric, non-finite, negative or huge value) and runs
# the backend on each. Every run must either succeed with an image or fail with
# the controlled validation exit code 2; a signal, a timeout, any other exit code,
# or an image from a failed run is reported.
# Usage: julia --project=. scripts/fuzz_native_requests.jl [max cases per request] [output directory]
using ROOTPlots

root_available() || error("local native ROOT is required")
limit = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 400
output = length(ARGS) >= 2 ? ARGS[2] : mktempdir()
mkpath(output)
config = ROOTPlots.find_root_config()

function run_case(exe, request, stem; timeout=60)
    process = run(pipeline(ignorestatus(Cmd([exe, request, stem])); stdout=devnull, stderr=devnull); wait=false)
    timed_out = Ref(false)
    timer = Timer(timeout) do _
        process_running(process) && (timed_out[] = true; kill(process, Base.SIGKILL))
    end
    wait(process); close(timer)
    png = isfile(stem * ".png")
    ok = !timed_out[] && process.termsignal == 0 &&
         ((process.exitcode == 0 && png) || (process.exitcode == 2 && !png))
    return (; ok, exitcode=process.exitcode, termsignal=process.termsignal, timed_out=timed_out[], png)
end

function mutations(valid)
    rows = split(chomp(valid), '\n')
    out = String[]
    for i in 2:length(rows)                                   # keep the schema header
        push!(out, join(deleteat!(copy(rows), i), '\n') * "\n")
        push!(out, join(insert!(copy(rows), i, rows[i]), '\n') * "\n")
        i < length(rows) && push!(out, join(rows[[1:i-1; i+1; i; i+2:end]], '\n') * "\n")
        fields = split(rows[i], '\t')
        for k in 2:length(fields), bad in ("", "x", "nan", "inf", "-1", "1e308", "99999999999")
            f = copy(fields); f[k] = bad
            push!(out, join([rows[1:i-1]; join(f, '\t'); rows[i+1:end]], '\n') * "\n")
        end
    end
    return unique(filter(!=(valid), out))
end

data = Hist1D([0.0, 1, 2, 3], [4.0, 6, 5], [4.0, 6, 5])
cmp = comparison(data, [Component("A", Hist1D([0.0, 1, 2, 3], [2.0, 3, 2], [2.0, 3, 2]), "#cc3311")])
g = GraphSeries([1.0, 2.0, 3.0], [2.0, 3.0, 2.5]; yerrors=[0.1, 0.2, 0.1], label="data")
requests = [
    ("comparison", ROOTPlots.build_root_backend(config),
     path -> ROOTPlots.write_root_spec(path, [cmp], [PanelSpec("x"; style=:CMS)], (0.5, 1.5))),
    ("cms graph", ROOTPlots.graphfit_executable(config),
     path -> (spec = GraphSpec(style=:CMS, cms_label=(extra_space=0.05, y_title_offset=1.1));
              ROOTPlots.write_graph_spec(path, [g], spec, ROOTPlots.graph_limits([g], spec)))),
]
failures = 0
for (name, exe, write_request) in requests
    path = joinpath(output, replace(name, ' ' => '_') * ".tsv")
    write_request(path)
    cases = mutations(read(path, String))
    n = min(limit, length(cases))
    for (k, case) in enumerate(cases[round.(Int, range(1, length(cases); length=n))])
        request = joinpath(output, "case.tsv")
        write(request, case)
        stem = joinpath(output, "out-$k")
        r = run_case(exe, request, stem)
        if !r.ok
            global failures += 1
            keep = joinpath(output, "failure-$(replace(name, ' ' => '_'))-$k.tsv")
            cp(request, keep; force=true)
            println("FAIL $name case $k: ", r, " (request kept: ", keep, ")")
        end
        rm(stem * ".png"; force=true); rm(stem * ".pdf"; force=true)
    end
    println(name, ": ", n, " of ", length(cases), " mutated requests run")
end
println(failures == 0 ? "no failures" : "$failures failures")
exit(failures == 0 ? 0 : 1)
