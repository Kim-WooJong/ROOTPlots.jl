using ROOTPlots: HistStack, Profile1D, Profile2D, Profile3D, PolygonHist2D, profile_stats

@testset "advanced histogram owned inputs and rejection" begin
    x=[.2,.4]; y=[1.,3.]; e=[0.,1.,2.]
    p=Profile1D(e,x,y)
    x[1]=.9; y[1]=99.; e[1]=-1.
    @test p.x==[.2,.4] && p.y==[1.,3.] && p.xedges==[0.,1.,2.]
    @test_throws ArgumentError Profile1D([0.,1.],[.2],[NaN])
    @test_throws ArgumentError Profile1D([0.,1.],[.2],[1.];weights=[-1.])
    @test_throws ArgumentError Profile1D([0.,1.],[.2],[1.];weights=[0.],error=:gaussian)
    @test_throws ArgumentError Profile1D([0.,1.],[.2],[1e308];weights=[2.])
    @test_throws ArgumentError Profile1D([0.,1.],[.2],[1.];error=:unknown)
    @test_throws ArgumentError Profile2D([0.,1.],[0.,1.],[.2],[.2,.3],[1.])
    @test_throws ArgumentError Profile3D([0.,1.],[0.,0.],[0.,1.],[.2],[.2],[.2],[1.])
    h=Hist1D([0.,1.,2.],[2.,3.],[2.,3.])
    s=HistStack([h,h];labels=["A","B"],fill_colors=["#cc3311","#0077bb"])
    h.sumw[1]=9.
    @test s.hists[1].sumw==[2.,3.]
    @test_throws ArgumentError HistStack([Hist1D([0.,1.],[-1.],[1.])])
    @test_throws ArgumentError HistStack([h,Hist1D([0.,2.,3.],[1.,1.],[1.,1.])])
    @test_throws ArgumentError HistStack([h];fill_styles=[-1])
    signed=HistStack([Hist1D([0.,1.],[-1.],[1.])];draw="NOSTACK HIST")
    @test signed.draw=="NOSTACK HIST"
    triangle=[0. 0.;1. 0.;.5 1.]
    poly=PolygonHist2D([triangle],[2.])
    triangle[1,1]=100.
    @test poly.polygons[1][1,1]==0.
    @test_throws ArgumentError PolygonHist2D([[0. 0.;1. 1.;0. 1.;1. 0.]],[1.])
    @test_throws ArgumentError PolygonHist2D([[0. 0.;1. 0.;2. 0.]],[1.])
    # Simple non-convex polygons are accepted; a repeated closing vertex is dropped.
    @test PolygonHist2D([[0. 0.;1. 0.;.2 .2;0. 1.]],[1.]).polygons[1]==[0. 0.;1. 0.;.2 .2;0. 1.]
    @test PolygonHist2D([[0. 0.;1. 0.;0. 1.;0. 0.]],[1.]).polygons[1]==[0. 0.;1. 0.;0. 1.]
    @test_throws ArgumentError PolygonHist2D([[0. 0.;1. 0.;1. 0.;0. 1.]],[1.])
    mktempdir() do dir
        missing=joinpath(dir,"not-created","plot")
        @test_throws ArgumentError render_root(signed,GraphSpec(),missing;draw="HIST",root_config="/not-root")
        @test_throws ArgumentError render_root(p,GraphSpec(logx=true),missing;root_config="/not-root")
        @test_throws ArgumentError render_root(p,GraphSpec(),missing;draw="E1 SAME",root_config="/not-root")
        @test_throws ArgumentError render_root(poly,Graph2DSpec(logx=true),missing;root_config="/not-root")
        @test_throws ArgumentError ROOTPlots.write_advanced_plot(joinpath(dir,"invalid.tsv"),p,Graph2DSpec())
        @test !isfile(joinpath(dir,"invalid.tsv"))
        @test !isdir(dirname(missing))
    end
end

@testset "native profile numerical definitions and multidimensional orientation" begin
    if TEST_NATIVE_ROOT
        p=Profile1D([0.,1.,2.],[.2,.4],[1.,3.])
        stats=profile_stats(p)
        @test stats.mean ≈ [2.,0.]
        @test stats.error ≈ [sqrt(.5),0.]
        @test stats.entries ≈ [2.,0.]
        @test stats.effective_entries ≈ [2.,0.]
        weighted=profile_stats(Profile1D([0.,1.],[.2,.4],[1.,3.];weights=[1.,3.]))
        @test weighted.mean[1]≈2.5
        @test weighted.error[1]≈sqrt(.75/1.6)
        @test weighted.entries[1]≈4.
        @test weighted.effective_entries[1]≈1.6
        spread=profile_stats(Profile1D([0.,1.],[.2,.4],[1.,3.];error=:spread))
        @test spread.error[1]≈1.
        gaussian=profile_stats(Profile1D([0.,1.],[.2,.4],[1.,3.];weights=[1.,3.],error=:gaussian))
        @test gaussian.error[1]≈.5
        integer=profile_stats(Profile1D([0.,1.],[.2],[2.];error=:integer))
        @test integer.error[1]≈1/sqrt(12.)
        @test profile_stats(Profile1D([0.,1.],Float64[],Float64[])).mean==[0.]
        # Lower edge belongs to the bin; upper edge and underflow are not visible.
        flow=profile_stats(Profile1D([0.,1.],[0.,1.,-.1],[2.,100.,200.]))
        @test flow.mean==[2.] && flow.entries==[1.]
        p2=Profile2D([0.,1.,2.],[0.,1.,2.],[.2,1.2],[1.2,.2],[3.,7.])
        stats2=profile_stats(p2)
        @test size(stats2.mean)==(2,2)
        @test stats2.mean==[0. 3.;7. 0.]
        p3=Profile3D([0.,1.,2.],[0.,1.,2.],[0.,1.,2.],[1.2],[.2],[1.2],[9.])
        stats3=profile_stats(p3)
        @test size(stats3.mean)==(2,2,2) && stats3.mean[2,1,2]==9.
        @test sum(stats3.entries)==1.
        overflow=Profile1D([0.,1.],[.2,.4],[1e154,1e154];error=:gaussian)
        redirect_stdio(stdout=devnull,stderr=devnull) do
            @test_throws ProcessFailedException profile_stats(overflow)
        end
    end
end

@testset "native advanced histogram graphics and bounded protocol" begin
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            p1=Profile1D([0.,1.,2.,3.],[.2,.4,1.2,1.4,2.2,2.4],[1.,3.,4.,5.,2.,3.])
            p2=Profile2D([0.,1.,2.],[0.,1.,2.],[.2,.4,1.2,1.4],[.2,1.2,.2,1.2],[1.,3.,4.,2.])
            p3=Profile3D([0.,1.,2.],[0.,1.,2.],[0.,1.,2.],
                [.2,1.2],[.2,1.2],[.2,1.2],[2.,5.])
            h=Hist1D([0.,1.,2.,3.],[2.,3.,1.],[2.,3.,1.])
            stack=HistStack([h,h];labels=["Signal","Background"],fill_colors=["#cc3311","#0077bb"])
            poly=PolygonHist2D([[0. 0.;1. 0.;.5 1.],[1. 0.;2. 0.;2. 1.;1. 1.]],[0.,4.])
            plots=((stack,GraphSpec(title="Stack")),(p1,GraphSpec(title="Profile")),
                (p2,Graph2DSpec()),(p2,Graph3DSpec()),(p3,Graph3DSpec(draw="BOX2Z")),(poly,Graph2DSpec()))
            for (i,(data,spec)) in enumerate(plots)
                paths=render_root(data,spec,joinpath(dir,"plot-$i");format=:both)
                @test filesize(paths.png)>2000
                @test read(paths.png,8)==UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a]
                @test read(paths.pdf,4)==collect(codeunits("%PDF"))
            end
            native=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            input=joinpath(dir,"bad.tsv")
            ROOTPlots.write_advanced_plot(input,p1,GraphSpec())
            original=read(input,String)
            for (i,broken) in enumerate((replace(original,"SAMPLE\t1\t"=>"SAMPLE\t99\t"),
                replace(original,"PEDGE\t1\t2\t1.0"=>"PEDGE\t1\t2\t0.0"),
                replace(original,"END\t6"=>"END\t9"),original*"EXTRA\n"))
                write(input,broken)
                process=run(pipeline(ignorestatus(Cmd([native,input,joinpath(dir,"bad-$i")]));stdout=devnull,stderr=devnull))
                @test process.exitcode==2
            end
            ROOTPlots.write_advanced_plot(input,poly,Graph2DSpec())
            write(input,replace(read(input,String),"VERTEX\t3\t0.5\t1.0"=>"VERTEX\t3\t0.5\t0.0"))
            process=run(pipeline(ignorestatus(Cmd([native,input,joinpath(dir,"bad-polygon")]));stdout=devnull,stderr=devnull))
            @test process.exitcode==2
            # Distinct-output process isolation remains compatible with rootmap.
            results=rootmap(1:2;workers=2) do i
                render_root(p1,GraphSpec(),joinpath(dir,"parallel-$i"))
            end
            @test all(r->isfile(r.png),results)
        end
    end
end
