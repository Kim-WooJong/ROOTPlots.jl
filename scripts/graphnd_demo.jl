using ROOTPlots

root_available() || error("local native ROOT is required for this demo")
format = isempty(ARGS) ? :png : Symbol(only(ARGS)) # png, pdf, or both
outdir = joinpath(@__DIR__, "..", "examples", "output")

# Synthetic scalar field sampled on rectangular cells (values[ix,iy]).
xedges = collect(range(-3.0, 3.0; length=25))
yedges = collect(range(-2.5, 2.5; length=21))
values = [4exp(-((xedges[ix]+xedges[ix+1])/2)^2/3 -
                ((yedges[iy]+yedges[iy+1])/2)^2/2) +
          0.5sin(2(xedges[ix]+xedges[ix+1])/2) for
          ix in 1:length(xedges)-1, iy in 1:length(yedges)-1]
grid = Grid2D(xedges, yedges, values)

println(render_root(grid, Graph2DSpec(title="Synthetic scalar field", xlabel="x", ylabel="y",
                                      zlabel="value"), joinpath(outdir,"graphnd_heatmap"); format))
println(render_root(grid, Graph2DSpec(title="Synthetic contours", xlabel="x", ylabel="y",
                                      draw=:contour, contours=12), joinpath(outdir,"graphnd_contour"); format))
println(render_root(grid, Graph3DSpec(title="Synthetic surface", xlabel="x", ylabel="y",
                                      zlabel="value"), joinpath(outdir,"graphnd_surface"); format))

t = collect(range(0, 4pi; length=90))
cloud = Graph3D(cos.(t), sin.(t), t ./ (2pi);
                style=GraphStyle(color="#c02731",marker_style=20,marker_size=0.8))
println(render_root(cloud, Graph3DSpec(title="Synthetic XYZ points", xlabel="x", ylabel="y",
                                       zlabel="z", draw=:points), joinpath(outdir,"graphnd_xyz"); format))
