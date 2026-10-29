@testset "special ROOT graphs and native primitives" begin
    @test_throws ArgumentError PolarGraph([0.,1.],[1.,-1.])
    @test_throws ArgumentError PolarGraph([0.,1.],[1.,2.];angle_errors=[0.,-1.])
    @test_throws ArgumentError QQGraph([1.],[1.])
    @test_throws ArgumentError SplineCurve([0.,0.,1.],[0.,1.,2.])
    @test_throws ArgumentError SplineCurve([0.,1.,2.],[0.,1.,2.];order=5)
    @test_throws ArgumentError PieChart([0.,0.])
    @test_throws ArgumentError PieChart([1.,2.];labels=["one"])
    @test_throws ArgumentError EfficiencyGraph([0.,1.],[2.],[1.])
    @test_throws ArgumentError EfficiencyGraph([0.,1.],[.5],[1.])
    @test_throws ArgumentError EfficiencyGraph([0.,1.],[1.],[2.];confidence=1.)
    source=[0.,1.,2.]; spline=SplineCurve(source,2 .* source;npx=3); source[1]=99.
    @test spline.x==[0.,1.,2.]
    g=GraphSeries([0.,1.,2.],[1.,2.,3.])
    primitives=PlotPrimitive[PlotArrow(.2,1.2,.8,1.8),PlotBox(.3,1.3,.5,1.5),
        PlotEllipse(1.,2.,.1,.2),SecondaryAxis(2.,1.,2.,3.,0.,100.;title="Percent")]
    for p in primitives
        @test startswith(ROOTPlots.primitive_record(p,GraphSpec()),"PRIMITIVE\t")
    end
    @test_throws ArgumentError ROOTPlots.primitive_record(PlotBox(1.,1.,0.,2.),GraphSpec())
    @test_throws ArgumentError ROOTPlots.primitive_record(PlotArrow(1.,1.,2.,2.;option="bad"),GraphSpec())
    @test_throws ArgumentError ROOTPlots.primitive_record(PlotEllipse(1.,1.,-1.,1.),GraphSpec())
    @test_throws ArgumentError ROOTPlots.primitive_record(SecondaryAxis(1.,1.,2.,2.,2.,1.),GraphSpec())
    @test_throws ArgumentError ROOTPlots.primitive_record(PlotBox(0.,1.,1.,2.),GraphSpec(logx=true))
    if TEST_NATIVE_ROOT
        efficiency=EfficiencyGraph([0.,1.,2.],[5.,10.],[10.,10.])
        points=efficiency_points(efficiency)
        quantiles=qq_points(QQGraph([1.,2.,3.,4.,5.],[2.,3.,6.,7.,10.]))
        @test quantiles.x==[1.,2.,3.,4.,5.]
        @test quantiles.y==[2.,3.,6.,7.,10.]
        @test points.x==[.5,1.5]
        @test points.y≈[.5,1.]
        @test points.eylow[1]≈points.eyhigh[1]
        @test 0<points.eylow[1]<.5 && points.eyhigh[2]==0
        @test spline_points(spline).y≈[0.,2.,4.]
        @test spline_points(SplineCurve(collect(0.:4.),collect(0.:4.);order=5,npx=5)).y≈collect(0.:4.)
        mktempdir() do dir
            polar=PolarGraph(collect(range(0.,2pi;length=30)),fill(1.,30))
            qq=QQGraph([1.,2.,3.,4.],[2.,4.,6.,8.])
            pie=PieChart([2.,3.,5.];labels=["A","B","C"],colors=["#4477aa","#cc6677","#228833"])
            bad=joinpath(dir,"absent","plot")
            @test_throws ArgumentError render_root(polar,GraphSpec(logx=true),bad)
            @test !isdir(dirname(bad))
            objects=[polar,qq,spline,pie,efficiency]
            files=rootmap(eachindex(objects);workers=3) do i
                render_root(objects[i],GraphSpec(title="Special $i"),joinpath(dir,"special-$i");format=:both)
            end
            @test all(p->read(p.png,8)==UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a],files)
            @test all(p->startswith(read(p.pdf,String),"%PDF-"),files)
            @test all(p->filesize(p.png)>1000,files)
            decorated=PlotPanel(g,GraphSpec(title="Primitives");primitives)
            @test isfile(render_root(decorated,joinpath(dir,"primitives")).png)
            panelgrid=[PlotPanel(object,GraphSpec(title="Special $i")) for (i,object) in enumerate(objects)]
            push!(panelgrid,decorated)
            @test isfile(render_root(PlotFigure(panelgrid;layout=(2,3)),joinpath(dir,"layout")).png)
            changed=EfficiencyGraph([0.,1.],[0.],[1.]); changed.passed[1]=2.
            @test_throws ArgumentError render_root(changed,GraphSpec(),bad)
        end
    end
end
