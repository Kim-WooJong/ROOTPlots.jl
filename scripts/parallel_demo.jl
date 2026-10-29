# julia --threads=4 --project=. scripts/parallel_demo.jl [output/directory]
# Synthetic independent measurements, not detector data.
using ROOTPlots

output = isempty(ARGS) ? joinpath(@__DIR__, "..", "examples", "output", "parallel") : ARGS[1]
x = collect(0.:0.25:5.)
jobs = [(id=i, data=GraphSeries(x, i .+ (i+1).*x;
                              yerrors=fill(0.2, length(x)), label="Synthetic data")) for i in 1:8]
build_root() # optional warm-up before the batch
results = rootmap(jobs; workers=4) do job
    fit = fit_root(job.data, :pol1; initial=[0., 0.])
    fit.status == 0 && fit.valid || error("fit failed for job $(job.id)")
    files = render_root([fit_curve(fit), job.data],
        GraphSpec(xlabel="x", ylabel="Response", title="Synthetic job $(job.id)"),
        joinpath(output, "job-$(job.id)"); format=:both)
    (; id=job.id, parameters=fit.parameters, files)
end
println("Julia default threads: ", Threads.nthreads(:default))
foreach(println, results)
