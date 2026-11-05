# Experiment styles: ROOT's built-in "ATLAS" TStyle, the package's "CMS" TStyle
# (cms-cat/cmsstyle values, cmsCanvas geometry and CMS_lumi label), and the CMS/ATLAS labels drawn with ROOT
# TLatex. The CMS/ATLAS labels are for collaboration members, under each
# collaboration's publication rules; the package only provides the styling.
using ROOTPlots

root_available() || error("local native ROOT is required for this demo")

# Reproducible pseudo-random samples (xorshift64 uniforms, Box-Muller Gaussians).
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

output = joinpath(@__DIR__, "..", "examples", "output", "experiment_styles")
files = String[]

# One histogram and one native ROOT Gaussian fit, drawn three times.
edges = collect(range(80.0, 100.0; length=41))
signal = histogram(91.2 .+ 2.5 .* gaussians(4000), edges).hist
fit = fit_root(signal, :gaus; initial=[400.0, 91.0, 2.5])
fit.status == 0 && fit.valid || error("ROOT fit did not converge")
common = (xlabel="m_{#mu#mu} [GeV]", ylabel="Events / 0.5 GeV", ylimits=(0.0, 520.0))

# 1. ROOT's default style (unchanged package look).
push!(files, render_root(signal, GraphSpec(; common...), joinpath(output, "fit_default");
                         draw=:e1, label="Data", attached_fit=fit, fit_label="Gaussian fit").png)
# 2. One switch: CMS style with the default "CMS Preliminary" label. As with cmsstyle's
# cmsCanvas on a non-square canvas, three-digit y labels need a larger y title offset
# than the 0.78 default and extra left space (cmsCanvas yTitOffset / extraSpace).
push!(files, render_root(signal, GraphSpec(; common..., style=:CMS,
                                           cms_label=(lumi="138 fb^{-1}", energy="13 TeV",
                                                      y_title_offset=1.1, extra_space=0.05)),
                         joinpath(output, "fit_cms"); draw=:e1, label="Data",
                         attached_fit=fit, fit_label="Gaussian fit").png)
# 3. ROOT's built-in ATLAS style with an "ATLAS Internal" label.
push!(files, render_root(signal, GraphSpec(; common..., style=:ATLAS,
                                           atlas_label=(text="Internal", energy="13 TeV", lumi="140 fb^{-1}")),
                         joinpath(output, "fit_atlas"); draw=:e1, label="Data",
                         attached_fit=fit, fit_label="Gaussian fit").png)
# 4. CMS fit/pull view: the label goes to the upper pad.
push!(files, render_root(signal, fit, GraphSpec(; common..., style=:CMS, legend=:topright,
                                                cms_label=(lumi="138 fb^{-1}", energy="13 TeV")),
                         joinpath(output, "fit_pull_cms"); data_label="Data").png)

# 5. CMS data vs. stacked prediction with a ratio panel.
mass = collect(range(60.0, 120.0; length=31))
bkg = histogram(60 .+ 60 .* uniforms(6000; seed=0x2545f4914f6cdd1d), mass).hist
sig = histogram(91.2 .+ 3.0 .* gaussians(9000; seed=0x0123456789abcdef), mass).hist
data = histogram(vcat(60 .+ 60 .* uniforms(6000; seed=0x1f2e3d4c5b6a7988),
                      91.2 .+ 3.0 .* gaussians(9000; seed=0x0fedcba987654321)), mass).hist
cmp = comparison(data, [Component("Background", bkg, "#6699cc"), Component("Z #rightarrow #mu#mu", sig, "#f2a33a")])
push!(files, render_root([cmp], [PanelSpec("m_{#mu#mu} [GeV]"; ylabel="Events / 2 GeV",
                                           lumi="138 fb^{-1} (13 TeV)", style=:CMS)],
                         joinpath(output, "stack_ratio_cms")).png)

# 6. CMS 2-D histogram (COLZ): the palette keeps its room on the right.
xy = histogram((gaussians(40000; seed=0x0123456789abcdef), 0.6 .* gaussians(40000; seed=0x0fedcba987654321)),
               (collect(range(-3.0, 3.0; length=41)), collect(range(-2.0, 2.0; length=31))))
push!(files, render_root(xy.hist, Graph2DSpec(xlabel="#eta", ylabel="#phi", zlabel="Events", style=:CMS,
                                              cms_label=(lumi="138 fb^{-1}", energy="13 TeV")),
                         joinpath(output, "hist2d_cms"); draw="COLZ").png)

# 7. CMS simulation label (no luminosity), placed outside the frame.
push!(files, render_root(signal, GraphSpec(; common..., style=:CMS,
                                           cms_label=(extra="Simulation", energy="13 TeV", position=:outside,
                                                      y_title_offset=1.1, extra_space=0.05)),
                         joinpath(output, "simulation_cms"); draw=:hist, label="Simulated signal").png)

# 8. 2x2 figure in CMS style: every pad gets its own label, sized for its pad.
figure = PlotFigure([
    PlotPanel(signal, GraphSpec(; common..., legend=:topright); draw=:e1, label="Data",
              attached_fit=fit, fit_label="Fit"),
    PlotPanel(xy.hist, Graph2DSpec(xlabel="#eta", ylabel="#phi", zlabel="Events"); draw="COLZ"),
    PlotPanel(bkg, GraphSpec(xlabel="m_{#mu#mu} [GeV]", ylabel="Events", legend=:none,
                             cms_label=(extra="Simulation",)); draw=:hist),
    PlotPanel([GraphSeries(collect(1.0:10.0), 2 .+ 0.3 .* collect(1.0:10.0); yerrors=fill(0.2, 10), label="Graph")],
              GraphSpec(xlabel="x", ylabel="y", legend=:bottomright, cms_label=false)),
]; layout=(2, 2), size=(1200, 900), style=:CMS)
push!(files, render_root(figure, joinpath(output, "figure_2x2_cms")).png)

# 9. CMS data vs. a four-process stacked prediction, log y, ratio panel 0.5-1.5:
# the ratio pad's text has the upper pad's size and "1.5"/"0.5" stay inside the pad.
pt = collect(35.0:10.0:205.0)
yields = (zjets=[90, 55, 20, 13, 12, 6.5, 5.5, 3.5, 2.5, 1.2, 2.0, 1.0, 0.6, 0.7, 0.6, 0.6, 2.5],
          diboson=[100, 95, 75, 62, 48, 38.5, 34.5, 26.5, 22.5, 18.8, 13, 12, 9.4, 7.3, 5.4, 3.4, 22.5],
          single_t=[210, 150, 155, 125, 110, 105, 85, 80, 75, 65, 60, 52, 45, 42, 39, 36, 35],
          ttbar=[2400, 1900, 1550, 1200, 930, 750, 575, 410, 280, 205, 145, 105, 85, 60, 40, 20, 200])
mc(v) = Hist1D(pt, Float64.(v), 0.02 .* v)   # MC with 2% statistical variance per bin
observed = [2600, 2100, 1750, 1350, 1050, 760, 600, 400, 300, 230, 170, 155, 120, 68, 65, 48, 220]
cmp = comparison(Hist1D(pt, Float64.(observed), Float64.(observed)),
                 [Component("Z + jets", mc(yields.zjets), "#2b7bba"), Component("Diboson", mc(yields.diboson), "#9467bd"),
                  Component("Single t", mc(yields.single_t), "#e377c2"), Component("t#bar{t}", mc(yields.ttbar), "#c00020")];
                 uncertainty_label="MC stat. unc.")
push!(files, render_root([cmp], [PanelSpec("Lepton p_{T} (GeV)"; ylabel="Events / 10 GeV", logy=true,
                                           ratio_limits=(0.5, 1.5), style=:CMS,
                                           cms_label=(extra="tutorial", lumi="1.21 fb^{-1}", energy="13.6 TeV",
                                                      position=:outside))],
                         joinpath(output, "stack_ratio_cms_logy")).png)

# 10. cmsstyle cmsCanvas(square=false, iPos=0, extraSpace=0.05, yTitOffset=1.1): markers
# without error bars, a dashed reference line in the legend, a right-hand text, log y.
charge = collect(0.5:1.0:12.5)
gain = [4.1e4, 7.9e4, 4.5e4, 7.2e4, 6.3e4, 8.7e4, 7.3e4, 6.3e4, 5.7e4, 7.0e4, 5.9e4, 4.6e4, 5.8e4]
mean_gain = sum(gain) / length(gain)
push!(files, render_root([GraphSeries(charge, gain; label="Measurement", style=GraphStyle(marker_style=20, marker_size=1.1))],
                         GraphSpec(xlabel="Integrated charge [C/cm^{2}]", ylabel="Effective gain",
                                   xlimits=(0.0, 13.0), ylimits=(5e3, 1e6), logy=true,
                                   lines=[PlotLine(0.0, mean_gain, 13.0, mean_gain; width=2, style=2,
                                                   label="Mean = $(round(mean_gain / 1e4; digits=2)) #times 10^{4}")],
                                   legend_spec=LegendSpec(bounds=(0.62, 0.75, 0.94, 0.88), text_size=0.04),
                                   style=:CMS, cms_label=(position=0, extra="p", right="University of Seoul",
                                                          extra_space=0.05, y_title_offset=1.1)),
                         joinpath(output, "stability_cms_outside")).png)

foreach(println, files)
