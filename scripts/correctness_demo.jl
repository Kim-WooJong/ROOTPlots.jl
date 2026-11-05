# Correctness fixes, drawn and reported by native ROOT:
#   1. log-y fit/pull views with empty tails and single-entry bins (native TH1D),
#      and log-y graphs whose error bars reach zero (clipped by ROOT);
#   2. Unicode labels translated to ROOT TLatex (unsupported text fails clearly);
#   3. native failures as RootBackendError with ROOT's own message;
#   4. ROOT Entries for histograms with unbinned moments but unknown fill count.
# Usage: julia --project=. scripts/correctness_demo.jl [output directory]
using ROOTPlots

root_available() || error("local native ROOT is required for this demo")
output = isempty(ARGS) ? joinpath(@__DIR__, "..", "examples", "output", "correctness") : ARGS[1]
mkpath(output)
files = String[]

# Reproducible Gaussian sample without a random-number package (xorshift64 + Box-Muller).
function gaussians(n; seed::UInt64=0x9e3779b97f4a7c15)
    state, u = seed, Vector{Float64}(undef, 2n)
    for i in eachindex(u)
        state ⊻= state << 13; state ⊻= state >> 7; state ⊻= state << 17
        u[i] = ((state >> 11) + 0.5) / 2.0^53
    end
    return [sqrt(-2log(u[2i-1])) * cospi(2u[2i]) for i in 1:n]
end

# 1. Log-y fit/pull: empty tail bins and bins with one entry (error bar down to zero).
h = histogram(gaussians(300), collect(range(-5.0, 5.0; length=51))).hist
println("empty bins: ", count(==(0), h.sumw), ", single-entry bins: ", count(==(1), h.sumw))
fit = fit_root(h, :gaus; initial=[20.0, 0.0, 1.0])
spec = GraphSpec(xlabel="x", ylabel="Entries / 0.2", logy=true)
push!(files, render_root(h, fit, spec, joinpath(output, "fit_pull_logy")).png)
push!(files, render_root(h, fit, GraphSpec(xlabel="x", ylabel="Entries / 0.2"), joinpath(output, "fit_pull_linear")).png)
g = GraphSeries([1.0, 2.0, 3.0, 4.0], [0.8, 12.0, 150.0, 2.0]; yerrors=[1.5, 4.0, 30.0, 2.0], label="errors reach zero")
push!(files, render_root([g], GraphSpec(xlabel="point", ylabel="value", logy=true, legend=:topleft),
                         joinpath(output, "graph_logy_clipped")).png)

# 2. Unicode labels -> TLatex.
for label in ("σ_m [GeV/c²] → Δφ", "α β", "ℓ⁺ℓ⁻ + E̸_T", "x̄ ± σ/√N, 5 °C, 10⁻³ cm⁻²")
    println(repr(label), " -> ", repr(ROOTPlots.root_latex(label)))
end
greek = GraphSpec(xlabel="σ_m [GeV/c²] → Δφ", ylabel="α β", title="Unicode labels: Γ(Z⁰ → ℓ⁺ℓ⁻)",
                  canvas=CanvasStyle(opt_title=1),
                  texts=[PlotText(0.15, 0.80, "αβγδεζηθικλμνξοπρστυφχψω"),
                         PlotText(0.15, 0.72, "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ"),
                         PlotText(0.15, 0.64, "± × ÷ · ° ≤ ≥ ≠ ≈ ∝ ∞ ∂ ∇ √ ∑ ∫ ∈ ⊂ ∩ ∪ ← → ⇒ ↔"),
                         PlotText(0.15, 0.56, "x² y₁₂ c⁻¹ x̄ ŷ p⃗ é ü Å ℏ ℓ")],
                  legend=:bottomright)
line = GraphSeries([0.0, 1.0, 2.0], [0.2, 0.25, 0.3]; label="μ⁺μ⁻ data", yerrors=[0.01, 0.01, 0.01])
push!(files, render_root([line], greek, joinpath(output, "unicode_labels")).png)
# The official cmsCanvas default on 800x600 (y title offset 0.78, no extra space)
# suits short y titles only: "Events / GeV" would reach the left canvas edge.
# cmsstyle's recipe for longer titles is cmsCanvas(extraSpace=0.05, yTitOffset=1.1),
# i.e. cms_label=(extra_space=0.05, y_title_offset=1.1); the CMS default is unchanged.
cms = GraphSpec(xlabel="m_{ℓℓ} [GeV]", ylabel="Events / GeV", style=:CMS,
                cms_label=(lumi="138 fb⁻¹", energy="13 TeV", additional_info="ℓ = e, μ",
                           extra_space=0.05, y_title_offset=1.1))
push!(files, render_root(h, cms, joinpath(output, "unicode_cms"); draw=:e1, label="μ events",
                         stats=true, stats_name="Δφ selection", stats_position=:topright).png)
try
    render_root([line], GraphSpec(xlabel="\uc9c8\ub7c9 [GeV]"), joinpath(output, "korean"))
catch err
    println("unsupported text: ", sprint(showerror, err))
end

# 3. Native failures keep ROOT's message, exit code and operation (also inside rootmap).
overflow = Profile1D([0.0, 1.0], [0.2, 0.4], [1e154, 1e154]; error=:gaussian)
try
    profile_stats(overflow)
catch err
    println(sprint(showerror, err))
end
try
    rootmap(p -> profile_stats(p).mean, [Profile1D([0.0, 1.0], [0.5], [1.0]), overflow]; workers=2)
catch err
    println(sprint(showerror, err))
end

# 4. Entries: moments without a fill count show ROOT's effective entries.
x = [-2.0, 0.25, 0.5, 0.75, 1.5, 1.75, 2.5, 3.0, 9.0]
w = [1.0, 2.0, 1.0, 0.5, 1.0, 1.5, 0.25, 1.0, 3.0]
filled = histogram(x, [0.0, 1.0, 2.0, 3.0]; weights=w).hist
unknown = Hist1D(filled.edges, filled.sumw, filled.sumw2; underflow=filled.underflow,
                 overflow=filled.overflow, moments=filled.moments)
fields = (:entries, :mean, :std, :underflow, :overflow)
push!(files, render_root(filled, GraphSpec(xlabel="x", ylabel="Σw"), joinpath(output, "entries_known");
                         draw=:e1, stats=fields, stats_name="entries known").png)
push!(files, render_root(unknown, GraphSpec(xlabel="x", ylabel="Σw"), joinpath(output, "entries_unknown");
                         draw=:e1, stats=fields, stats_name="moments, entries unknown").png)

println("wrote:")
foreach(f -> println("  ", f), files)
