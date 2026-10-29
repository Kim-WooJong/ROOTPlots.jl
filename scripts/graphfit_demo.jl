using ROOTPlots

root_available() || error("local native ROOT is required for this demo")
x = collect(0.0:4.0)
y = [1.1, 2.9, 5.2, 6.8, 9.1] # synthetic observations
observed = GraphSeries(x, y; yerrors=fill(0.2, length(x)),
                       label="Synthetic measurements")
fit = fit_root(observed, :pol1; initial=[0.0, 1.0])
fit.status == 0 && fit.valid || error("ROOT fit did not converge")
files = render_root([fit_curve(fit), observed],
                    GraphSpec(xlabel="x", ylabel="Response", title="ROOT graph and fit",
                              legend=:topleft),
                    joinpath(@__DIR__, "..", "examples", "output", "graphfit_demo"))
println("parameters = ", fit.parameters)
println("errors = ", fit.errors)
println("covariance = ", fit.covariance)
println("status = ", fit.status, ", chi2/ndf = ", fit.chi2, "/", fit.ndf)
println(files)
