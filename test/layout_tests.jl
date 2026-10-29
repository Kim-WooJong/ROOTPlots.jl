@testset "native multi-pad layout contracts and ND errors" begin
    h=Hist1D([0.,1.,2.],[2.,4.],[2.,4.])
    g=GraphSeries([0.,1.,2.],[1.,3.,2.];yerrors=[.1,.2,.1],label="data")
    h2=Hist2D([0.,1.,2.],[0.,1.,2.],[1. 3.;2. 4.],fill(.25,2,2))
    h3=Hist3D([0.,1.,2.],[0.,1.,2.],[0.,1.,2.],reshape(collect(1.:8.),2,2,2),fill(.5,2,2,2))
    panel=PlotPanel([g],GraphSpec(title="Graph",legend=:topleft))
    @test_throws ArgumentError PlotFigure(PlotPanel[])
    @test_throws ArgumentError PlotFigure([panel];layout=(0,1))
    @test_throws ArgumentError PlotFigure([panel,panel];layout=(1,1))
    @test_throws ArgumentError PlotFigure([panel];size=(0,100))
    matrix=[panel PlotPanel(h,GraphSpec(title="Hist"));panel panel]
    @test PlotFigure(matrix).layout==(2,2)
    @test PlotFigure(matrix).panels[2].spec.title=="Hist"
    mktempdir() do dir
        bad=joinpath(dir,"absent","bad")
        @test_throws ArgumentError render_root(PlotFigure([PlotPanel(h,GraphSpec();unknown=true)]),bad)
        @test_throws ArgumentError render_root(PlotFigure([PlotPanel(h2,Graph2DSpec();contour_levels=[1.,1.])]),bad)
        @test_throws ArgumentError render_root(PlotFigure([PlotPanel(h2,Graph2DSpec();zlimits=(2.,1.))]),bad)
        @test_throws ArgumentError render_root(PlotFigure([
            PlotPanel(h2,Graph2DSpec(palette=:bird)),PlotPanel(h2,Graph2DSpec(palette=:viridis))]),bad)
        @test_throws ArgumentError render_root(PlotFigure([
            PlotPanel(h,GraphSpec();canvas=CanvasStyle(end_error_size=1.)),
            PlotPanel(h2,Graph2DSpec();canvas=CanvasStyle(end_error_size=2.))]),bad)
        @test !isdir(dirname(bad))
        payload=ROOTPlots.write_panel_request(joinpath(dir,"panel.tsv"),PlotPanel(h2,Graph2DSpec(contours=17);
            zaxis=AxisStyle(title_size=.06,center_title=true),zlimits=(.1,5.),contour_levels=[1.,2.,3.]))
        text=read(payload,String)
        @test occursin("CELL\t1\t1\t1.0\t0.25",text)
        @test occursin("LEVELS\t3\t1.0\t2.0\t3.0",text)
        legend_spec=GraphSpec(legend_spec=LegendSpec(entries=[LegendEntry(:series,2,"second")]))
        multi=ROOTPlots.write_panel_request(joinpath(dir,"multi.tsv"),PlotPanel([g,g],legend_spec;
            labels=["first","overridden"],styles=[GraphStyle(),GraphStyle(color="#cc3311")]))
        @test occursin("overridden\t#cc3311",read(multi,String))
        if TEST_NATIVE_ROOT
            stack=HistStack([h,h];labels=["A","B"],fill_colors=["#cc3311","#0077bb"])
            profile=Profile1D([0.,1.,2.],[.2,.4,1.2,1.4],[1.,3.,2.,5.])
            curve=ModelCurve(:gaus,[1.,0.,1.];range=(-3.,3.))
            xyz=Graph3D([0.,1.,0.,1.],[0.,0.,1.,1.],[1.,2.,3.,4.])
            poly=PolygonHist2D([([0.,1.,.5],[0.,0.,1.])],[2.])
            panels=[panel,PlotPanel(h,GraphSpec(title="Histogram")),
                PlotPanel(h2,Graph2DSpec(title="TH2");zaxis=AxisStyle(center_title=true)),
                PlotPanel(h3,Graph3DSpec(title="TH3",draw="BOX2Z")),
                PlotPanel(stack,GraphSpec(title="Stack")),PlotPanel(profile,GraphSpec(title="Profile")),
                PlotPanel(curve,GraphSpec(title="TF1")),PlotPanel(xyz,Graph3DSpec(title="XYZ",draw=:points)),
                PlotPanel(poly,Graph2DSpec(title="Polygon"))]
            files=render_root(PlotFigure(panels;layout=(3,3),size=(1200,1000)),joinpath(dir,"all");format=:both)
            @test read(files.png,8)==UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a]
            @test startswith(read(files.pdf,String),"%PDF-")
            @test filesize(files.png)>10000
            @test isfile(render_root(h2,Graph2DSpec(),joinpath(dir,"h2");draw="COLZ",zlimits=(0.,5.)).png)
            @test isfile(render_root(h3,Graph3DSpec(draw="BOX2Z"),joinpath(dir,"h3")).png)
            # A trailing XY pad must not change the previous ND pad's palette.
            bird=PlotPanel(h2,Graph2DSpec(palette=:bird))
            xy=PlotPanel(h,GraphSpec())
            @test isfile(render_root(PlotFigure([bird,xy]),joinpath(dir,"bird-then-xy")).png)
            styled=render_root(h,GraphSpec(),joinpath(dir,"filled");draw="HIST",
                style=GraphStyle(fill_color="#cc3311",fill_alpha=.8))
            @test isfile(styled.png)
            @test isfile(render_root([h,h],GraphSpec(),joinpath(dir,"filled-overlay");
                draw=["HIST","HIST SAME"],styles=[GraphStyle(fill_color="#cc3311"),GraphStyle(fill_color="#4477aa")]).png)
            @test isfile(render_root(PlotPanel([g,ModelCurve(:pol1,[1.,1.];range=(0.,2.))],GraphSpec()),
                joinpath(dir,"actual-function-overlay")).png)
            exe=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            malformed=joinpath(dir,"malformed.tsv")
            write(malformed,"ROOTPlotsGraph\t1\nDRAWLAYOUT\t800\t600\t1\t1\t1\nPANEL\t1\t999\nEND\t1\n")
            process=run(pipeline(ignorestatus(`$exe $malformed $(joinpath(dir,"bad-native"))`);stdout=devnull,stderr=devnull))
            @test process.exitcode==2
            copies=rootmap(1:2;workers=2) do i
                render_root(PlotFigure(panels[1:2]),joinpath(dir,"parallel-$i"))
            end
            @test all(p->isfile(p.png),copies)
        end
    end
end
