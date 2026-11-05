# Safe histogram operations, checked in ROOT's own statistics boxes:
#   1. two independent samples, each scaled (cross section x luminosity / events),
#      then merged - identical boxes to one direct fill with the scaled weights;
#   2. 2-D -> 1-D projections: default (visible bins only), flow=true with exact
#      statistics, and flow=true where integrated-axis flows force ROOT's fallback;
#   3. a 3-D -> 2-D projection with flow=true.
# The event values are also written as text so a plain ROOT macro can fill
# TH1D/TH2D with them (TH1::Scale, TH1::Add, TH2::ProjectionX) for comparison.
# Usage: julia --project=. scripts/histogram_ops_demo.jl [output directory]
using ROOTPlots

root_available() || error("local native ROOT is required for this demo")
output = isempty(ARGS) ? joinpath(@__DIR__, "..", "examples", "output", "histogram_ops") : ARGS[1]
mkpath(output)
files = String[]

# Reproducible Gaussian samples without a random-number package (xorshift64 + Box-Muller).
function gaussians(n; seed::UInt64)
    state, u = seed, Vector{Float64}(undef, 2n)
    for i in eachindex(u)
        state ⊻= state << 13; state ⊻= state >> 7; state ⊻= state << 17
        u[i] = ((state >> 11) + 0.5) / 2.0^53
    end
    return [sqrt(-2log(u[2i-1])) * cospi(2u[2i]) for i in 1:n]
end

"""Render with ROOTPLOTS_STATS_DUMP and print the box lines ROOT drew."""
function render_with_box(f, name)
    dump = joinpath(output, name * ".stats.txt")
    rm(dump; force=true)
    path = withenv(f, "ROOTPLOTS_STATS_DUMP" => dump)
    lines = [split(l, '\t')[2] for l in readlines(dump) if startswith(l, "LINE\t")]
    println(rpad(name, 22), join(lines, " | "))
    push!(files, path)
end

fields = (:name, :entries, :mean, :std, :underflow, :overflow)
edges = collect(range(-4.0, 4.0; length=41))
xa = 0.5 .+ 1.0 .* gaussians(3000; seed=0x9e3779b97f4a7c15)
xb = -0.3 .+ 1.4 .* gaussians(1500; seed=0x2545f4914f6cdd1d)
wa, wb = 0.8, 1.5                          # per-sample scale factors
open(joinpath(output, "samples_1d.txt"), "w") do io
    foreach(x -> println(io, "a\t", x), xa); foreach(x -> println(io, "b\t", x), xb)
end

# 1. Scale each independent sample, then merge; compare with one direct weighted fill.
ha, hb = histogram(xa, edges).hist, histogram(xb, edges).hist
merged = merge(wa * ha, wb * hb)
direct = histogram([xa; xb], edges; weights=[fill(wa, length(xa)); fill(wb, length(xb))]).hist
spec = GraphSpec(xlabel="x", ylabel="Weighted events / 0.2", legend=:none)
println("statistics boxes drawn by ROOT:")
render_with_box("merged_scaled") do
    render_root(merged, spec, joinpath(output, "merged_scaled"); draw=:hist, stats=fields,
                stats_name="0.8 A + 1.5 B (merge)").png
end
render_with_box("direct_fill") do
    render_root(direct, spec, joinpath(output, "direct_fill"); draw=:hist, stats=fields,
                stats_name="direct weighted fill").png
end
println("max |bin difference| merged vs direct: ", maximum(abs.(merged.sumw .- direct.sumw)))
try
    merge(ha, ha)
catch err
    println("merge(ha, ha): ", sprint(showerror, err))
end

# 2. 2-D -> 1-D projections. y is correlated with x; the narrow y range has flows.
gx, gy = gaussians(5000; seed=0x0123456789abcdef), gaussians(5000; seed=0x7f4a7c159e3779b9)
x2, y2 = gx, 0.6 .* gx .+ 0.8 .* gy
open(joinpath(output, "samples_2d.txt"), "w") do io
    foreach(i -> println(io, x2[i], '\t', y2[i]), eachindex(x2))
end
e = collect(range(-2.5, 2.5; length=21))
wide = collect(range(-6.0, 6.0; length=25))
narrow2 = histogram((x2, y2), (e, e)).hist          # y flows under visible x: fallback
wide2 = histogram((x2, y2), (e, wide)).hist          # no y flows: exact statistics
pspec = GraphSpec(xlabel="x", ylabel="Events / 0.25", legend=:none)
for (name, h, flow, title) in (("projx_default", narrow2, false, "default (visible bins)"),
                               ("projx_flow_exact", wide2, true, "flow=true, exact"),
                               ("projx_flow_fallback", narrow2, true, "flow=true, y flows"))
    p = projection(h, :x; flow)
    println(rpad(name, 22), "entries=", p.entries, " moments ", p.moments === nothing ? "unknown" : "known")
    render_with_box(name) do
        render_root(p, pspec, joinpath(output, name); draw=:hist, stats=fields, stats_name=title).png
    end
end
direct_x = histogram(x2, e).hist
render_with_box("x_direct_fill") do
    render_root(direct_x, pspec, joinpath(output, "x_direct_fill"); draw=:hist, stats=fields,
                stats_name="direct fill of x").png
end

# 3. 3-D -> 2-D (z integrated with its flows); the result keeps x/y flow bins.
z3 = gaussians(5000; seed=0x5851f42d4c957f2d)
h3 = histogram((x2, y2, z3), (e, wide, collect(range(-5.0, 5.0; length=11)))).hist
pxy = projection(h3, (:x, :y); flow=true)
println("3D -> (x, y) flow=true: entries=", pxy.entries, ", moments known: ", pxy.moments !== nothing)
render_with_box("proj3d_xy") do
    render_root(pxy, Graph2DSpec(draw=:heatmap, xlabel="x", ylabel="y", size=(700, 600)),
                joinpath(output, "proj3d_xy"); stats=true).png
end

println("wrote:")
foreach(f -> println("  ", f), files)
