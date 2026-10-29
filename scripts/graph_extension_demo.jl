using ROOTPlots

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
