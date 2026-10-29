using ROOTPlots

# Synthetic inputs only. Writes PNG and PDF files to examples/output/graph-extensions/.
outdir=joinpath(@__DIR__,"..","examples","output","graph-extensions")
mkpath(outdir)
x=collect(0.5:0.25:4.)
y=2 .+ sin.(x)
band=GraphSeries(x,y;yerrors=(low=fill(.2,length(x)),high=fill(.35,length(x))),
                 errors=:band,draw=:line,label="Asymmetric errors",
                 style=GraphStyle(color="#336699",fill_alpha=.35))
curve=function_curve(t->2+sin(t),(.5,4.);n=200,label="Julia callable samples",
                     style=GraphStyle(color="#cc3311"))
println(render_root([band,curve],GraphSpec(xlabel="x",ylabel="response"),
                    joinpath(outdir,"band-and-callable");format=:both))
println(render_root(ModelCurve(:gaus,[3.,2.,.5];range=(0.,4.),label="Native TF1"),
                    GraphSpec(xlabel="x",ylabel="Gaussian"),joinpath(outdir,"native-tf1");format=:both))
xx=[0.,1.,0.,1.,.3,.7]; yy=[0.,0.,1.,1.,.6,.4]; zz=1 .+ xx.^2 .+ yy
cloud=Graph3D(xx,yy,zz)
println(render_root(cloud,Graph3DSpec(draw=:surface,interpolate=true,resolution=(40,40),
                    title="Explicit ROOT scattered interpolation"),joinpath(outdir,"scattered-surface");format=:both))
errors=Graph3DErrors(xx,yy,zz;zerrors=fill(.15,length(xx)))
println(render_root(errors,Graph3DSpec(draw=:points,title="Native XYZ error bars"),
                    joinpath(outdir,"xyz-errors");format=:both))

# TGraphBentErrors: two series measured at the same x; bending the y-error bars
# left/right keeps both visible without moving the points.
px=collect(1.:6.)
a=BentErrorSeries(px,[2.,2.6,3.1,3.0,3.6,4.1];yerrors=(low=fill(.45,6),high=fill(.35,6)),
                  ybend=fill(-.18,6),label="Series A (bent left)",style=GraphStyle(color="#cc3311",marker_style=21))
b=BentErrorSeries(px,[2.1,2.5,3.0,3.2,3.5,4.0];yerrors=fill(.4,6),
                  ybend=fill(.18,6),label="Series B (bent right)",style=GraphStyle(color="#336699",marker_style=20))
println(render_root([a,b],GraphSpec(xlabel="x",ylabel="y",title="TGraphBentErrors",legend=:topleft),
                    joinpath(outdir,"bent-errors");format=:both))

# TGraphMultiErrors: statistical bars plus systematic boxes, styled per set.
stat=YErrorSet(fill(.25,6);label="Statistical",style=GraphStyle(color="#000000"))
syst=YErrorSet((low=[.3,.35,.5,.45,.4,.6],high=[.4,.3,.45,.55,.5,.5]);label="Systematic",
               display=:boxes,xscale=.6,style=GraphStyle(color="#4477aa",fill_alpha=.45))
multi=MultiErrorSeries(px,[2.,2.6,3.1,3.0,3.6,4.1],[stat,syst];xerrors=fill(.3,6),label="Synthetic data")
println(render_root(multi,GraphSpec(xlabel="x",ylabel="y",title="TGraphMultiErrors",legend=:topleft),
                    joinpath(outdir,"multi-errors");format=:both))
println("combined (quadrature) errors: ",total_yerrors(multi))

# Native TF2 from a predefined formula: surface, contours, and heatmap pads.
surface=ModelSurface(:bigaus,[1.,0.,1.,0.,.7,.5];xrange=(-3.,3.),yrange=(-3.,3.),npx=100,npy=100)
mesh=ModelSurface(:bigaus,[1.,0.,1.,0.,.7,.5];xrange=(-3.,3.),yrange=(-3.,3.),npx=30,npy=30)
println(render_root(PlotFigure([
        PlotPanel(mesh,Graph3DSpec(draw="SURF1",title="TF2 SURF1",zlabel="f(x,y)")),
        PlotPanel(surface,Graph2DSpec(draw="CONT1Z",contours=12,title="TF2 CONT1Z")),
        PlotPanel(surface,Graph2DSpec(draw=:heatmap,title="TF2 COLZ"))];layout=(1,3),size=(1500,500)),
    joinpath(outdir,"tf2");format=:both))

# Native TF3: implicit surface f = level, and natively sampled TH3D boxes.
volume=ModelVolume(:gaus3d,[1.,0.,1.,0.,.7,0.,1.3];xrange=(-3.,3.),yrange=(-3.,3.),zrange=(-3.,3.),level=.3)
println(render_root(PlotFigure([
        PlotPanel(volume,Graph3DSpec(draw=:iso,title="TF3 iso-surface f = 0.3")),
        PlotPanel(volume,Graph3DSpec(draw=:box,title="TF3 sampled into TH3D, bins with f > 0.3 (BOX2Z)"))];size=(1200,550)),
    joinpath(outdir,"tf3");format=:both))

# A Julia callable sampled on a grid (display-only TH2D, not a TF2).
ripple=function_grid((u,v)->sin(u)*cos(v)*exp(-(u^2+v^2)/20),(-6.,6.),(-6.,6.);n=(120,120))
println(render_root(ripple,Graph2DSpec(title="Julia callable on a grid",zlabel="f"),
                    joinpath(outdir,"function-grid");format=:both))

# Non-convex polygons: TH2Poly bins and filled polygon primitives.
L=([0.,2,2,1,1,0],[0.,0,1,1,3,3])
star=([3.2+(isodd(i) ? .45 : 1.1)*cos(pi/2+i*pi/5) for i in 0:9],
      [1.6+(isodd(i) ? .45 : 1.1)*sin(pi/2+i*pi/5) for i in 0:9])
U=([2.5,4.5,4.5,4.,4.,3.,3.,2.5],[3.,3.,4.5,4.5,3.5,3.5,4.5,4.5])
cells=PolygonHist2D([L,star,U],[1.,3.,2.])
frame=GraphSeries([0.,4.6],[0.,4.6];draw=:points,style=GraphStyle(marker_size=.01))
shapes=PlotPanel([frame],GraphSpec(title="PlotPolygon primitives",legend=:none);
                 primitives=[PlotPolygon(L...;fill_color="#cc6677"),PlotPolygon(star...;fill_color="#228833",width=2),
                             PlotPolygon(U...;fill_color="#4477aa",fill_alpha=.3)])
println(render_root(PlotFigure([PlotPanel(cells,Graph2DSpec(title="Non-convex TH2Poly bins")),shapes];size=(1200,550)),
                    joinpath(outdir,"nonconvex-polygons");format=:both))

# Nested pads: a 2x1 sub-grid inside the right-hand cell, with weighted widths.
right=PlotFigure([PlotPanel(surface,Graph2DSpec(draw=:heatmap,title="Sub-grid: TF2")),
                  PlotPanel([a,b],GraphSpec(title="Sub-grid: bent errors",legend=:none))];
                 layout=(2,1),heights=[1,1])
println(render_root(PlotFigure([PlotPanel(multi,GraphSpec(title="Outer cell: multi errors",legend=:topleft)),right];
                               layout=(1,2),widths=[3,2],size=(1300,700)),
                    joinpath(outdir,"nested-layout");format=:both))
