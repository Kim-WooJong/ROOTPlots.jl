# ROOT axis log scales, titles, and number/label display: AxisStyle, CanvasStyle
# and LabelChange fields map one-to-one onto ROOT setters (TAxis, TPad, gStyle).
# Each figure compares an option with ROOT's default where that is instructive.
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

output = joinpath(@__DIR__, "..", "examples", "output", "axis_format")
files = String[]
small = (600, 450)

# 1. Log x and log y on a 1-D histogram (TPad::SetLogx/SetLogy), log-spaced bins.
energy = 10 .^ (3 .* uniforms(20_000) .^ 2)                       # falling spectrum on 1..1000
spectrum = histogram(energy, 10 .^ range(0, 3; length=31)).hist
# A larger x title offset keeps the title clear of the 10^3 label at the axis end.
push!(files, render_root(spectrum, GraphSpec(xlabel="Energy (keV)", ylabel="Entries / bin",
    title="Log x and log y", logx=true, logy=true, legend=:none, xaxis=AxisStyle(title_offset=1.3)),
    joinpath(output, "01_log_xy"); draw="E1").png)

# 2. Log z on a 2-D COLZ histogram (TPad::SetLogz); the palette is logarithmic.
gx, gy = gaussians(50_000; seed=0x1111222233334444), gaussians(50_000; seed=0x5555666677778888)
xy = histogram((gx, 0.6 .* gx .+ 0.8 .* gy), (collect(-4.0:0.2:4.0), collect(-4.0:0.2:4.0))).hist
push!(files, render_root(xy, Graph2DSpec(title="Log z on COLZ", xlabel="x", ylabel="y", zlabel="Entries",
    logz=true), joinpath(output, "02_logz_colz"); draw="COLZ").png)

# 3. more_log_labels (TAxis::SetMoreLogLabels) vs ROOT's default on a short log range.
rate = GraphSeries([1., 2., 3., 4., 5., 6.], [2.2, 4.1, 7.9, 15.0, 31.0, 62.0]; label="rate")
logspec(title; yaxis=AxisStyle()) = GraphSpec(; xlabel="Step", ylabel="Rate (Hz)", title, logy=true,
    legend=:none, yaxis, ylimits=(1.5, 90.))
push!(files, render_root(PlotFigure([
    PlotPanel([rate], logspec("ROOT default")),
    PlotPanel([rate], logspec("more_log_labels=true"; yaxis=AxisStyle(more_log_labels=true)))];
    layout=(1, 2), size=(1000, 450)), joinpath(output, "03_more_log_labels")).png)

# 4. no_exponent (TAxis::SetNoExponent) vs ROOT's default x10^n with large numbers.
counts = GraphSeries(collect(1.:6.), [1.2e5, 2.6e5, 4.1e5, 3.3e5, 2.0e5, 0.9e5]; draw=:linepoints, label="counts")
bigspec(title; yaxis=AxisStyle()) = GraphSpec(; xlabel="Run", ylabel="Events", title, legend=:none, yaxis,
    canvas=CanvasStyle(left_margin=0.16))
push!(files, render_root(PlotFigure([
    PlotPanel([counts], bigspec("ROOT default")),
    PlotPanel([counts], bigspec("no_exponent=true"; yaxis=AxisStyle(no_exponent=true)))];
    layout=(1, 2), size=(1000, 450)), joinpath(output, "04_no_exponent")).png)

# 5. decimals (TAxis::SetDecimals) and divisions (TAxis::SetNdivisions(n, optim)).
eff = GraphSeries(collect(0.:0.5:3.), [0.2, 0.55, 0.8, 0.9, 1.0, 1.05, 1.1]; draw=:linepoints, label="eff")
fracspec(title; xaxis=AxisStyle(), yaxis=AxisStyle()) = GraphSpec(; xlabel="Threshold", ylabel="Efficiency",
    title, legend=:none, xaxis, yaxis, xlimits=(0., 3.), ylimits=(0., 1.2))
push!(files, render_root(PlotFigure([
    PlotPanel([eff], fracspec("ROOT default")),
    PlotPanel([eff], fracspec("decimals=true"; yaxis=AxisStyle(decimals=true))),
    PlotPanel([eff], fracspec("divisions=3 / 12, not optimized";
        xaxis=AxisStyle(divisions=3, optimize_divisions=false),
        yaxis=AxisStyle(divisions=12, optimize_divisions=false)))];
    layout=(1, 3), size=(1350, 450)), joinpath(output, "05_decimals_divisions")).png)

# 6. Bin labels (TAxis::SetBinLabel) drawn vertically (TAxis::LabelsOption("v")).
channels = ["e#mu", "ee", "#mu#mu", "e#tau", "#mu#tau", "#tau#tau", "jets", "other"]
yields = Hist1D(collect(0.:8.), [120., 80., 95., 40., 35., 12., 300., 18.], [120., 80., 95., 40., 35., 12., 300., 18.])
push!(files, render_root(yields, GraphSpec(xlabel="", ylabel="Events", title="Bin labels, vertical", legend=:none,
    xaxis=AxisStyle(bin_labels=channels, labels_option="v", label_size=0.05),
    canvas=CanvasStyle(bottom_margin=0.16)), joinpath(output, "06_bin_labels"); draw="HIST").png)

# 7. Time axis: x values are seconds after time_offset (TAxis::SetTimeDisplay/SetTimeFormat/SetTimeOffset).
seconds = collect(0.0:900.0:86_400.0)
temperature = GraphSeries(seconds, 18 .+ 4 .* sinpi.((seconds .- 21_600) ./ 43_200); draw=:line, label="T")
push!(files, render_root([temperature], GraphSpec(xlabel="Time (GMT), 2025-01-01", ylabel="Temperature (#circC)",
    title="Time axis", legend=:none,
    xaxis=AxisStyle(time_display=true, time_format="%H:%M", time_offset=1_735_689_600, divisions=508)),
    joinpath(output, "07_time_axis")).png)

# 8. change_labels (TAxis::ChangeLabel): recolour/rotate the first label, rename one, hide the last.
line = GraphSeries(collect(0.:10.), collect(0.:10.) .^ 1.5; draw=:line, label="f")
push!(files, render_root([line], GraphSpec(xlabel="x", ylabel="f(x)", title="change_labels", legend=:none,
    xaxis=AxisStyle(change_labels=[LabelChange(1; color="#cc3311", angle=30, size=0.05),
                                   LabelChange(3; text="four"), LabelChange(-1; size=0)])),
    joinpath(output, "08_change_labels")).png)

# 9. Grid lines (TPad::SetGridx/SetGridy) with ticks on both sides of the axes (TAxis::SetTicks("+-")).
push!(files, render_root([rate], GraphSpec(xlabel="Step", ylabel="Rate (Hz)", title="Grid lines and +- ticks",
    legend=:none, logy=true, xaxis=AxisStyle(ticks="+-", title_offset=1.3), canvas=CanvasStyle(grid_x=true, grid_y=true)),
    joinpath(output, "09_grid")).png)

# 10. Title: the package's TLatex heading (default) vs ROOT's own title box (gStyle OptTitle).
#     With opt_title=true, ROOT's "title;x title;y title" form also sets the axis titles.
gauss = histogram(gaussians(5000), collect(-4.0:0.25:4.0)).hist
push!(files, render_root(PlotFigure([
    PlotPanel(gauss, GraphSpec(title="Default heading", xlabel="x", ylabel="Entries", legend=:none); draw=:e1),
    PlotPanel(gauss, GraphSpec(title="ROOT title box;x (from title);Entries (from title)", legend=:none,
        canvas=CanvasStyle(opt_title=true)); draw=:e1),
    PlotPanel(gauss, GraphSpec(title="Styled ROOT title", xlabel="x", ylabel="Entries", legend=:none,
        canvas=CanvasStyle(opt_title=true, title_x=0.12, title_align=13, title_border_size=1,
                           title_fill_color="#f2f2f2", title_style=1001, title_text_color="#003366",
                           title_font_size=0.06)); draw=:e1)];
    layout=(1, 3), size=(1350, 450)), joinpath(output, "10_title_box")).png)

# 11. Axis colours (TAxis::SetAxisColor/SetLabelColor/SetTitleColor).
push!(files, render_root([counts], GraphSpec(xlabel="Run", ylabel="Events", title="Axis colours", legend=:none,
    xaxis=AxisStyle(axis_color="#0077bb", label_color="#0077bb", title_color="#0077bb"),
    yaxis=AxisStyle(axis_color="#cc3311", label_color="#cc3311", title_color="#cc3311", max_digits=3)),
    joinpath(output, "11_axis_colours")).png)

# 12. 2x2 figure with different axis settings in every pad.
push!(files, render_root(PlotFigure([
    PlotPanel(spectrum, GraphSpec(xlabel="Energy (keV)", ylabel="Entries", title="log x/y, more labels", legend=:none,
        logx=true, logy=true, xaxis=AxisStyle(more_log_labels=true, title_offset=1.3), canvas=CanvasStyle(grid_x=true)); draw=:e1),
    PlotPanel(yields, GraphSpec(xlabel="", ylabel="Events", title="Bin labels", legend=:none,
        xaxis=AxisStyle(bin_labels=channels, labels_option="u", label_size=0.06), canvas=CanvasStyle(bottom_margin=0.15)); draw=:hist),
    PlotPanel(xy, Graph2DSpec(title="log z, ROOT title", zlabel="Entries", logz=true,
        zaxis=AxisStyle(label_color="#0077bb"), canvas=CanvasStyle(opt_title=true)); draw="COLZ"),
    PlotPanel([temperature], GraphSpec(xlabel="Time (GMT)", ylabel="#circC", title="Time axis", legend=:none,
        xaxis=AxisStyle(time_display=true, time_format="%Hh", time_offset=1_735_689_600),
        yaxis=AxisStyle(decimals=true)))];
    layout=(2, 2), size=(1100, 850)), joinpath(output, "12_figure_2x2")).png)

# 13. 2-D histogram with a styled z (palette) axis and large counts.
big = Hist2D(xy.xedges, xy.yedges, 250 .* xy.sumw, 250^2 .* xy.sumw2)
push!(files, render_root(big, Graph2DSpec(title="z axis style", xlabel="x", ylabel="y", zlabel="Entries",
    zaxis=AxisStyle(no_exponent=true, divisions=505, label_color="#555555", title_offset=1.6, label_size=0.03),
    canvas=CanvasStyle(right_margin=0.22)), joinpath(output, "13_zaxis_style"); draw="COLZ").png)

# 14. Two-pad views: comparison ratio pad and fit/pull layout with axis styles and grids.
data = histogram(energy[1:2000], 10 .^ range(0, 3; length=11)).hist
model = Hist1D(data.edges, data.sumw .* 1.05, data.sumw .* 1.05)
cmp = comparison(data, [Component("Model", model, "#4477aa")])
push!(files, render_root([cmp], [PanelSpec("Energy (keV)"; logx=true, logy=true,
    xaxis=AxisStyle(more_log_labels=true), ratio_yaxis=AxisStyle(divisions=504, decimals=true),
    canvas=CanvasStyle(grid_y=true))], joinpath(output, "14_ratio_axes")).png)
fit = fit_root(gauss, :gaus; initial=[300., 0., 1.])
push!(files, render_root(gauss, fit, GraphSpec(xlabel="x", ylabel="Entries", title="Fit / pull",
    xaxis=AxisStyle(divisions=508), canvas=CanvasStyle(grid_x=true, opt_title=true)),
    joinpath(output, "15_fit_pull_axes")).png)

foreach(println, files)
