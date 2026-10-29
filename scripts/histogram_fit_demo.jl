using DataFrames, ROOTPlots

root_available() || error("local native ROOT is required for this demo")

# Synthetic events: DataFrames owns selection and derived columns.
counts = [8, 12, 15, 11]
mass = vcat([fill((i - 0.5) / 2, n) for (i, n) in enumerate(counts)]...)
events = DataFrame(mass=vcat(mass, 100.0), accepted=vcat(trues(length(mass)), false))
selected = subset(events, :accepted => ByRow(identity))
transform!(selected, :mass => ByRow(x -> 2x) => :observable)

# ROOTPlots owns binning, ROOT fitting, graphics, and fit diagnostics.
h = histogram(selected, :observable, collect(0.0:4.0)).hist
fit = fit_root(h, :pol0; initial=[sum(h.sumw) / length(h.sumw)])
fit.status == 0 && fit.valid || error("ROOT fit did not converge")

output = joinpath(@__DIR__, "..", "examples", "output")
files = render_root(h, fit,
                    GraphSpec(xlabel="Observable", ylabel="Events", legend=:topleft,
                              title="Histogram fit and pulls"),
                    joinpath(output, "histogram_fit_demo"); format=:both,
                    data_label="Binned events")

println("fit χ²/ndf = ", fit.chi2, "/", fit.ndf)
println("χ² p-value = ", chi2_pvalue(fit))
println("fit/pull plot = ", files)
