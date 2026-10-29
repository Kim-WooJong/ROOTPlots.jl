#!/usr/bin/env julia
# Run: julia --project=. scripts/advanced_histogram_demo.jl [output-directory]
using ROOTPlots

function advanced_histogram_demo(outdir=joinpath(@__DIR__,"..","plots","advanced"))
    mkpath(outdir)
    h1=Hist1D([0.,1.,2.,3.],[4.,7.,3.],[4.,7.,3.])
    h2=Hist1D([0.,1.,2.,3.],[2.,3.,5.],[2.,3.,5.])
    stack=HistStack([h1,h2];labels=["Signal","Background"],fill_colors=["#cc3311","#0077bb"])
    files=Dict{Symbol,Any}()
    files[:stack]=render_root(stack,GraphSpec(title="Native THStack",xlabel="x",ylabel="Counts"),
                              joinpath(outdir,"stack");format=:both)
    x=[.2,.4,1.2,1.4,2.2,2.4]
    p1=Profile1D([0.,1.,2.,3.],x,[1.,3.,4.,5.,2.,3.])
    files[:profile1d]=render_root(p1,GraphSpec(title="Native TProfile",ylabel="Mean response"),
                                  joinpath(outdir,"profile1d");format=:both)
    p2=Profile2D([0.,1.,2.],[0.,1.,2.],[.2,.4,1.2,1.4],[.2,1.2,.2,1.2],[1.,3.,4.,2.])
    files[:profile2d]=render_root(p2,Graph2DSpec(title="Native TProfile2D",zlabel="Mean response"),
                                  joinpath(outdir,"profile2d");format=:both)
    files[:profile2d_surface]=render_root(p2,Graph3DSpec(title="Profile surface"),
                                          joinpath(outdir,"profile2d-surface");format=:both)
    p3=Profile3D([0.,1.,2.],[0.,1.,2.],[0.,1.,2.],[.2,1.2],[.2,1.2],[.2,1.2],[2.,5.])
    files[:profile3d]=render_root(p3,Graph3DSpec(draw="BOX2Z",title="Native TProfile3D"),
                                  joinpath(outdir,"profile3d");format=:both)
    poly=PolygonHist2D([[0. 0.;1. 0.;.5 1.],[1. 0.;2. 0.;2. 1.;1. 1.]],[1.,4.])
    files[:polygon]=render_root(poly,Graph2DSpec(title="Native TH2Poly"),
                                joinpath(outdir,"polygon");format=:both)
    println("Native profile summary: ",profile_stats(p1))
    return files
end

abspath(PROGRAM_FILE)==@__FILE__ && advanced_histogram_demo(isempty(ARGS) ? joinpath(@__DIR__,"..","plots","advanced") : ARGS[1])
