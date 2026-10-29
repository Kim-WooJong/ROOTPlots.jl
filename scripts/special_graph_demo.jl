using ROOTPlots
# Synthetic examples: native ROOT polar, QQ, spline, pie, efficiency and shapes.
theta=collect(range(0.,2pi;length=80))
graphs=[PolarGraph(theta,1 .+ .3cos.(3theta)),
    QQGraph([1.,2.,3.,4.,5.],[2.,3.,6.,7.,10.]),
    SplineCurve(collect(0.:5.),[0.,1.,.5,2.,1.5,3.]),
    PieChart([2.,3.,5.];labels=["A","B","C"],colors=["#4477aa","#cc6677","#228833"]),
    EfficiencyGraph([0.,1.,2.,3.,4.],[1.,5.,8.,10.],fill(10.,4))]
names=["TGraphPolar","TGraphQQ","TSpline3","TPie","TEfficiency"]
panels=[PlotPanel(g,GraphSpec(title=names[i])) for (i,g) in enumerate(graphs)]
g=GraphSeries([0.,1.,2.],[1.,2.,3.])
push!(panels,PlotPanel(g,GraphSpec(title="ROOT annotation primitives");
    primitives=[PlotArrow(.2,2.5,1.,2.),PlotBox(.3,1.3,.6,1.6;fill_style=1001),
                PlotEllipse(1.,2.,.2,.3),SecondaryAxis(2.,1.,2.,3.,0.,100.;title="Percent")]))
println(render_root(PlotFigure(panels;layout=(2,3),size=(1400,900)),
    joinpath(@__DIR__,"..","examples","output","special-graphs");format=:both))
