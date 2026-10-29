# ROOT-style statistics and fit boxes (TPaveStats) drawn by native ROOT.
using ROOTPlots

root_available() || error("local native ROOT is required for this demo")

# Reproducible pseudo-random samples (xorshift64 uniforms, Box-Muller Gaussians);
# no random-number package is needed.
function uniforms(n; seed::UInt64=0x9e3779b97f4a7c15)
    state, out = seed, Vector{Float64}(undef, n)
    for i in 1:n
        state ⊻= state << 13; state ⊻= state >> 7; state ⊻= state << 17
        out[i] = ((state >> 11) + 0.5) / 2.0^53
    end
    return out
end
function gaussians(n; seed::UInt64=0x9e3779b97f4a7c15)
    u = uniforms(2n; seed)
    return [sqrt(-2log(u[2i-1])) * cospi(2u[2i]) for i in 1:n]
end
values = 0.3 .+ 1.1 .* gaussians(5000)
append!(values, [-6.0, -5.5, 7.0, 8.0, 9.5])          # a few under/overflow entries

output = joinpath(@__DIR__, "..", "examples", "output", "stats_box")
spec = GraphSpec(xlabel="Observable", ylabel="Entries")

# histogram() records Entries (every fill, flows included), flows, and unbinned moments,
# so ROOT shows the same Entries/Mean/Std Dev as a TH1D filled event by event.
h = histogram(values, collect(range(-4.0, 4.0; length=41))).hist
files = String[]

# 1. ROOT's default box: name, Entries, Mean, Std Dev.
push!(files, render_root(h, spec, joinpath(output, "hist_default"); draw="E1", label="signal",
                         stats=true).png)

# 2. Gaussian fit drawn as a native TF1 with ROOT's fit box (χ²/ndf, Prob, parameters).
#    The tall box would cover the peak; stats_headroom=true raises the y maximum
#    so the data stay below it (ROOT's default, false, draws the box over the data).
fit = fit_root(h, :gaus; initial=[200.0, 0.0, 1.0])
fit.status == 0 && fit.valid || error("ROOT fit did not converge")
push!(files, render_root(h, GraphSpec(xlabel="Observable", ylabel="Entries", legend=:topleft),
                         joinpath(output, "hist_fit"); draw=:e1, label="signal",
                         stats=true, fitstats=true, attached_fit=fit, fit_label="Gaussian fit",
                         stats_headroom=true).png)

# 3. Custom fields, including under/overflow and errors on mean/std dev (ROOT letters
#    "nemrou" would select the same fields without errors).
push!(files, render_root(h, GraphSpec(xlabel="Observable", ylabel="Entries", legend=:topleft),
                         joinpath(output, "hist_fields"); draw=:hist, label="signal",
                         stats=(:name, :entries, :mean_error, :std_error, :underflow,
                                :overflow, :integral, :skewness), stats_headroom=true).png)

# 4. 2-D histogram: Entries, Mean x/y, Std Dev x/y.
xy = histogram((gaussians(20000; seed=0x0123456789abcdef), 0.6 .* gaussians(20000; seed=0x0fedcba987654321)),
               (collect(range(-3.0, 3.0; length=31)), collect(range(-2.0, 2.0; length=21))))
push!(files, render_root(xy.hist, Graph2DSpec(xlabel="x", ylabel="y", zlabel="Entries"),
                         joinpath(output, "hist2d"); stats=true, label="xy").png)

# 5. Profile: ROOT's profile box (Mean, Mean y, Std Dev, Std Dev y).
px = 4 .* uniforms(2000; seed=0x5555aaaa5555aaaa)
py = 2 .* px .+ 0.5 .* gaussians(2000; seed=0x1234123412341234)
profile = Profile1D(collect(0.0:0.25:4.0), px, py)
push!(files, render_root(profile, GraphSpec(xlabel="x", ylabel="<y>", legend=:none),
                         joinpath(output, "profile"); stats=true, stats_name="response",
                         stats_position=:bottomright).png)

# 6. Graph with a fit: ROOT shows only the fit box for graphs.
gx = collect(1.0:1.0:10.0)
gy = 1.0 .+ 2.0 .* gx .+ 0.4 .* gaussians(10; seed=0x0badc0ffee0ddf00)
graph = GraphSeries(gx, gy; yerrors=fill(0.4, 10), label="Measurements")
gfit = fit_root(graph, :pol1; initial=[0.0, 1.0])
push!(files, render_root([graph], GraphSpec(xlabel="x", ylabel="y", legend=:topleft),
                         joinpath(output, "graph_fit"); attached_fit=gfit, fitstats=true,
                         fit_label="Linear fit", stats_position=:bottomright).png)

# 7. Boxes inside a 2x2 PlotFigure: each pad keeps its own box settings.
figure = PlotFigure([
    PlotPanel(h, GraphSpec(xlabel="Observable", ylabel="Entries", legend=:topleft); draw=:e1,
              label="signal", stats=true, attached_fit=fit, fitstats=true, stats_headroom=true),
    PlotPanel(xy.hist, Graph2DSpec(xlabel="x", ylabel="y"); label="xy", stats=true),
    PlotPanel(profile, GraphSpec(xlabel="x", ylabel="<y>", legend=:none); stats=true,
              stats_position=:bottomright),
    PlotPanel([graph], GraphSpec(xlabel="x", ylabel="y", legend=:topleft); attached_fit=gfit,
              fitstats=true, stats_position=:bottomright),
]; layout=(2, 2), size=(1200, 900))
push!(files, render_root(figure, joinpath(output, "figure_2x2")).png)

foreach(println, files)
