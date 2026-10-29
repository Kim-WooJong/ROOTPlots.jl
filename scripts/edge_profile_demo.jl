# Publication-style edge fits on synthetic data (native ROOT required).
#   julia --project=. scripts/edge_profile_demo.jl [output-dir]
# Data are generated deterministically; they are not measurements.
using ROOTPlots

root_available() || error("native ROOT is required; see docs/installation.md")
outdir = isempty(ARGS) ? joinpath(@__DIR__, "..", "examples", "output", "edge-profile") : ARGS[1]

"""Deterministic pseudo-noise in [-1, 1] (integer hash) so the example is reproducible."""
wiggle(i) = 2 * ((UInt32(i) * 0x9E3779B1 ⊻ 0x5bd1e995) * 0x85ebca6b % 0x00100000) / 0x00100000 - 1

publication_axis() = AxisStyle(title_font=42, label_font=42, title_size=0.055,
                               label_size=0.045, title_offset=1.1, divisions=505)
publication_canvas() = CanvasStyle(left_margin=0.15, right_margin=0.04, bottom_margin=0.14,
                                   top_margin=0.05, ticks_x=true, ticks_y=true,
                                   frame_line_width=2, end_error_size=0.)

# 1. A single falling edge (step response) fitted with erfc_edge = [B, A, R, sigma].
x = collect(range(0.05, 2.95; length=40))
truth = 0.2 .+ (5.0 / 2) .* [ROOTPlots.erfc((xi - 1.5) / (sqrt(2) * 0.18)) for xi in x]
errors = 0.05 .+ 0.03 .* sqrt.(truth)
y = truth .+ 0.6 .* errors .* wiggle.(eachindex(x))
edge = GraphSeries(x, y; yerrors=errors, label="Synthetic profile")
fit = fit_root(edge, :erfc_edge; initial=[0.1, 4.0, 1.4, 0.2],
               bounds=[nothing, nothing, nothing, (0.02, 1.0)])
fit.status == 0 && fit.valid || error("edge fit failed")
R, sigma = fit.parameters[3], fit.parameters[4]
spec = GraphSpec(xlabel="Position x (mm)", ylabel="Signal (a.u.)", size=(900, 720),
    xaxis=publication_axis(), yaxis=publication_axis(), canvas=publication_canvas(),
    lines=[PlotLine(R, 0., R, maximum(y); color="#666666", style=2, label="Fitted edge")],
    texts=[PlotText(0.60, 0.48, "R = $(round(R; digits=3)) #pm $(round(fit.errors[3]; digits=3))"; size=0.036),
           PlotText(0.60, 0.42, "#sigma = $(round(sigma; digits=3)) #pm $(round(fit.errors[4]; digits=3))"; size=0.036)])
edge_files = render_root([fit_curve(fit; label="erfc fit"), edge], spec,
                         joinpath(outdir, "edge"); format=:both)

# 2. A box profile with two edges fitted with double_erf = [B, H, xL, xR, sigmaL, sigmaR].
xb = collect(range(-3.0, 3.0; length=61))
box = [10 + 300 / 2 * (ROOTPlots.erf((v + 1.2) / (sqrt(2) * 0.25)) -
                       ROOTPlots.erf((v - 1.2) / (sqrt(2) * 0.30))) for v in xb]
h = Hist1D(collect(range(-3.05, 3.05; length=62)), round.(box .+ 2 .* sqrt.(box) .* wiggle.(eachindex(xb))),
           max.(box, 1.0))
boxfit = fit_root(h, :double_erf; initial=[10., 280., -1.0, 1.0, 0.3, 0.3],
                  bounds=[nothing, nothing, nothing, nothing, (0.02, 1.0), (0.02, 1.0)])
boxfit.status == 0 && boxfit.valid || error("box fit failed")
box_files = render_root(h, boxfit, GraphSpec(xlabel="Position x (mm)", ylabel="Counts / bin",
                                             legend=:topright),
                        joinpath(outdir, "box-fit-pull"); format=:both)

println("edge: R=", R, " sigma=", sigma, " chi2/ndf=", fit.chi2, "/", fit.ndf)
println("box:  xL=", boxfit.parameters[3], " xR=", boxfit.parameters[4],
        " chi2/ndf=", boxfit.chi2, "/", boxfit.ndf)
println(edge_files); println(box_files)
