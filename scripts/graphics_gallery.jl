# Synthetic ROOT graphics gallery. No ROOT-file I/O or RDataFrame.
# julia --threads=4 --project=. scripts/graphics_gallery.jl
using ROOTPlots
x=collect(range(0.1,5.;length=100))
data=GraphSeries(x,2 .+ sin.(x);yerrors=fill(.15,length(x)),label="Synthetic data",errors=:band,
                 style=GraphStyle(fill_color="#4477aa",fill_alpha=.3))
curve=function_curve(t->2+sin(t),(0.1,5.);label="Julia callable")
h=histogram([.2,.4,.6,1.2,1.4,1.8,2.2,2.4,2.8,2.9],[0.,1.,2.,3.]).hist
stack=HistStack([h,h];labels=["A","B"],fill_colors=["#cc6677","#4477aa"])
profile=Profile1D([0.,1.,2.,3.],[.2,.4,1.2,1.4,2.2,2.4],[1.,2.,3.,5.,4.,6.])
h2=Hist2D([0.,1.,2.],[0.,1.,2.],[1. 2.;3. 4.],fill(.25,2,2))
h3=Hist3D([0.,1.,2.],[0.,1.,2.],[0.,1.,2.],reshape(collect(1.:8.),2,2,2),ones(2,2,2))
poly=PolygonHist2D([([0.,1.,.5],[0.,0.,1.]),([1.,2.,1.5],[0.,0.,1.])],[2.,5.])
xyz=Graph3D([0.,1.,0.,1.],[0.,0.,1.,1.],[1.,2.,3.,5.])
panels=[PlotPanel([data,curve],GraphSpec(title="Graph + band + callable")),
    PlotPanel(h,GraphSpec(title="TH1 errors");draw=:e1),
    PlotPanel(stack,GraphSpec(title="THStack")),
    PlotPanel(profile,GraphSpec(title="TProfile")),
    PlotPanel(h2,Graph2DSpec(title="TH2 + stored variance");zaxis=AxisStyle(title_size=.06),zlimits=(0.,5.)),
    PlotPanel(h3,Graph3DSpec(title="TH3",draw="BOX2Z")),
    PlotPanel(poly,Graph2DSpec(title="TH2Poly")),
    PlotPanel(ModelCurve(:gaus,[1.,0.,1.];range=(-3.,3.)),GraphSpec(title="TF1")),
    PlotPanel(xyz,Graph3DSpec(title="TGraph2D",draw=:surface,interpolate=true))]
figure=PlotFigure(panels;layout=(3,3),size=(1500,1200))
println(render_root(figure,joinpath(@__DIR__,"..","examples","output","graphics-gallery");format=:both))
