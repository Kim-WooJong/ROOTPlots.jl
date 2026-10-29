include("setup.jl")
using DataFrames

@testset "ROOTPlots" begin
include("dependency_tests.jl")
include("scene_tests.jl")
include("build_tests.jl")
include("parallel_tests.jl")
include("histogramnd_tests.jl")
include("graph_extension_tests.jl")
include("advanced_histogram_tests.jl")
include("layout_tests.jl")
include("special_graph_tests.jl")
include("graph_feature_tests.jl")
include("visual_tests.jl")

@testset "DataFrame histogram: ROOT edges, signed weights, flow bins" begin
    df = DataFrame(x=[-1.0, 0.0, 1.0, 2.0], w=[2.0, 1.0, -3.0, 4.0])
    result = histogram(df, :x, [0.0, 1.0, 2.0]; weight=:w)
    @test result.hist.sumw == [1.0, -3.0]
    @test result.hist.sumw2 == [1.0, 9.0]
    @test result.underflow == (2.0, 4.0)
    @test result.overflow == (4.0, 16.0)
    @test_throws ArgumentError Hist1D([0, 0, 1], [1, 2], [1, 2])
end

@testset "explicit whole-bin rebinning" begin
    h = Hist1D([0.,1,3,6], [2.,-1,4], [2.,5,10])
    merged = rebin(h, [0.,3,6])
    @test merged.edges == [0.,3,6]
    @test merged.sumw == [1.,4]
    @test merged.sumw2 == [7.,10]
    @test integral(merged) == integral(h)
    @test rebin(h, h.edges).sumw == h.sumw
    merged.sumw[1] = 100.
    merged.edges[2] = 2.
    @test h.sumw == [2.,-1,4] && h.edges == [0.,1,3,6]
    for bad in ([0.,1.5,6], [0.,3], [1.,3,6], [0.,3,7],
                [0.,6,3], [0.,0,6])
        @test_throws ArgumentError rebin(h, bad)
    end
    damaged = Hist1D([0.,1,2],[1.,2.],[1.,2.])
    damaged.sumw2[2] = -1.
    @test_throws ArgumentError rebin(damaged, [0.,2.])

    data = Hist1D([0.,1,3,6], [2.,1,4], [2.,1,4])
    mc = Hist1D([0.,2,3,6], [1.,2,3], [1.,2,3])
    @test_throws ArgumentError comparison(data, [Component("MC", mc, "#336699")])
    cmp = comparison(rebin(data, [0.,3,6]),
                     [Component("MC", rebin(mc, [0.,3,6]), "#336699")])
    @test cmp.prediction == [3.,3]
    @test cmp.ratio == [1.,4/3]
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            files = render_root([cmp], [PanelSpec("x")], joinpath(dir,"rebinned"))
            @test isfile(files.png) && filesize(files.png) > 1000
        end
    end
end

@testset "visible-bin statistics and fit diagnostics" begin
    h = Hist1D([0.,1,2,3,4], [1.,2,3,4], [1.,2,3,4])
    @test integral(h).value == 10
    @test integral(h).error ≈ sqrt(10)
    @test integral(h; bins=2:3) == (value=5.0, error=sqrt(5.0))
    @test effective_entries(h) ≈ 10
    @test binned_mean(h) == 2.5
    @test binned_std(h) == 1
    @test binned_mean_error(h) ≈ 1/sqrt(10)
    s = histstats(h; underflow=(2.,4.), overflow=(3.,9.))
    @test s.mean == 2.5 && s.integral.value == 10
    @test s.underflow == (2.,4.) && s.overflow == (3.,9.)
    @test !s.negative_bins
    variable = Hist1D([0.,1,2,4], [2.,3,1], [2.,3,1])
    @test binned_mean(variable) ≈ 17/12
    @test binned_std(variable) ≈ sqrt(101)/12
    @test binned_mean_error(variable) ≈ 0.3419037118080881 # ROOT 6.40 TH1::GetMeanError
    @test integral(variable; width=true).value == 7
    @test integral(variable; width=true).error ≈ 3
    points = histogram_points(variable; label="Binned data")
    @test points.x == [0.5, 1.5, 3.0]
    @test points.y == [2.0, 3.0, 1.0]
    @test points.eylow == points.eyhigh == sqrt.([2.0, 3.0, 1.0])
    @test all(iszero, points.exlow) && all(iszero, points.exhigh)
    points.y[1] = 99.0
    @test variable.sumw[1] == 2.0
    steps = histogram_steps(variable; label="Binned line")
    @test steps.x == [0.,1,1,2,2,4]
    @test steps.y == [2.,2,3,3,1,1]
    @test steps.draw == :line && steps.errors == :none && steps.kind == :plain
    steps.y[1] = 99.0
    @test variable.sumw[1] == 2.0
    @test histogram_steps(Hist1D([0.,1.],[-2.],[4.])).y == [-2.,-2.]
    invalid_steps = Hist1D([0.,1.],[1.],[1.])
    invalid_steps.edges[2] = 0.
    @test_throws ArgumentError histogram_steps(invalid_steps)
    @test_throws ArgumentError integral(h; bins=[1,1])
    @test_throws ArgumentError histstats(h; underflow=(1.,-1.))
    @test_throws ArgumentError binned_mean(Hist1D([0.,1.],[0.],[0.]))
    @test effective_entries(Hist1D([0.,1.],[0.],[0.])) == 0
    signed = Hist1D([0.,1,2], [3.,-1.], [5.,1.])
    @test binned_mean(signed) == 0
    @test_throws ArgumentError binned_std(signed)
    @test_throws ArgumentError histstats(signed)
    @test_throws ArgumentError binned_mean(Hist1D([0.,1,2], [1.,-1.], [1.,1.]))

    mktempdir() do dir
        path = joinpath(dir,"result.tsv")
        write(path,"ROOTPlotsFitResult\t1\nFIT\t0\t1\t3\t3.841458820694124\t1\t2\t0.0\tMinuit2\t0\n" *
                   "PARAM\t1\t1.0\t0.5\t0\nCOV\t1\t0.25\nEND\n")
        f = ROOTPlots.read_fit_result(path,:pol0,(0.,2.),:chi2)
        @test chi2_pvalue(f) ≈ 0.05 atol=1e-10
        @test fit_correlation(f) == ones(1,1)
        r = fit_residuals(Hist1D([0.,1,2],[2.,0.],[2.,0.]), f)
        @test r.valid == [true,false]
        @test r.residual[1] == 1 && r.pull[1] ≈ 1/sqrt(2)
        @test isnan(r.pull[2])
        pull_graph = fit_pull_points(Hist1D([0.,1,2],[2.,0.],[2.,0.]), f)
        @test pull_graph.x == [0.5] && pull_graph.y ≈ [1/sqrt(2)]
        @test_throws ArgumentError fit_pull_points(Hist1D([0.,1,2],[0.,0.],[0.,0.]), f)
        mock_h = Hist1D([0.,1,2],[2.,0.],[2.,0.])
        badstem = joinpath(dir,"not-created","plot")
        @test_throws ArgumentError render_root(mock_h,f,GraphSpec(),badstem;
                                               pull_limits=(0.,2.))
        @test_throws ArgumentError render_root(mock_h,f,GraphSpec(size=(10001,600)),badstem)
        @test_throws ArgumentError render_root(mock_h,f,GraphSpec(),badstem;format=:bad)
        @test !isdir(dirname(badstem))
        weighted = ROOTPlots.read_fit_result(path,:pol0,(0.,2.),:weighted_poisson)
        @test weighted.chi2_kind == :weighted_poisson_objective
        @test_throws ArgumentError chi2_pvalue(weighted)
        @test_throws ArgumentError render_root(mock_h,weighted,GraphSpec(),badstem)
    end
end

@testset "ROOT graph input and axis validation" begin
    s = GraphSeries([2,1,3], [4,2,6]; yerrors=(low=[0.2,0.3,0.4],high=[0.5,0.6,0.7]),
                    label="Observed",draw=:linepoints)
    @test s.kind == :asymmetric
    @test s.x == [2,1,3] # connected graphs preserve caller order
    @test ROOTPlots.graph_limits([s],GraphSpec())[1][1] < 1
    @test_throws ArgumentError GraphSeries([1,2],[3])
    @test_throws ArgumentError GraphSeries([1],[2];yerrors=[-1])
    @test_throws ArgumentError GraphSeries([1],[2];label="bad\tlabel")
    @test_throws ArgumentError GraphSeries([1],[2];style=GraphStyle(color="red"))
    @test_throws ArgumentError ROOTPlots.graph_limits([GraphSeries([0,1],[1,2])],
                                                       GraphSpec(logx=true,xlimits=(0.1,2.0)))
    @test_throws ArgumentError ROOTPlots.graph_limits([GraphSeries([1,2],[0,2])],
                                                       GraphSpec(logy=true,ylimits=(0.1,3.0)))
    @test_throws ArgumentError ROOTPlots.graph_limits([GraphSeries([1],[2];xerrors=[1])],
                                                       GraphSpec(logx=true))
    @test_throws ArgumentError ROOTPlots.graph_limits([GraphSeries([1e308,1.79e308],[1,2])],
                                                       GraphSpec())
end

@testset "native TH1 drawing input and protocol" begin
    h = Hist1D([0.,1,2,4], [2.,3,1], [2.,3,1])
    spec = GraphSpec(xlabel="x", ylabel="Count", legend=:none)
    @test ROOTPlots.hist_draw_mode("HIST BAR") == :bar
    @test ROOTPlots.hist_draw_mode("e1") == :e1
    @test ROOTPlots.hist_draw_mode("P") == :points
    @test ROOTPlots.hist_draw_mode("TEXT") == :text
    @test ROOTPlots.hist_draw_mode(:steps) == :steps
    @test ROOTPlots.hist_draw_limits(h,spec,:e1,GraphStyle(),"Histogram")[2][2] >
          ROOTPlots.hist_draw_limits(h,spec,:hist,GraphStyle(),"Histogram")[2][2]
    badstem = joinpath(tempdir(),"rootplots-invalid-hist-option","plot")
    for bad in (:unknown, "HIST SAME", "HIST\nE1", "E1\tX0", repeat("A",128),
                "é", "SPEC BF(1600", "SPEC \";evil()", "SPEC \\bad", "HIST;do_something()",
                "GLBOX")
        @test_throws ArgumentError render_root(h,spec,badstem;draw=bad)
    end
    @test ROOTPlots.root_draw_option("Pie") == "Pie"
    @test ROOTPlots.root_draw_option("") == ""
    @test_throws ArgumentError render_root(h,spec,badstem;draw="FUNC")
    @test !isdir(dirname(badstem))
    @test_throws ArgumentError render_root(h,GraphSpec(logx=true),badstem;draw=:e1)
    @test_throws ArgumentError render_root(h,GraphSpec(logy=true,ylimits=(-1.,10.)),
                                           badstem;draw=:e2)
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            exe = ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            input = joinpath(dir,"hist.tsv")
            ROOTPlots.write_hist_draw_spec(input,h,spec,
                ROOTPlots.hist_draw_limits(h,spec,:hist,GraphStyle(),"H"),
                :hist,"H",GraphStyle())
            valid = read(input,String)
            for malformed in (
                replace(valid,"\thist\n"=>"\tunsupported\n";count=1),
                replace(valid,"EDGE\t2\t1.0"=>"EDGE\t2\t0.0";count=1),
                replace(valid,"BIN\t1\t2.0\t2.0"=>"BIN\t1\t2.0\t-1.0";count=1),
                valid * "END\t1\n",
            )
                write(input,malformed)
                @test !success(pipeline(Cmd([exe,input,joinpath(dir,"bad")]);
                                        stdout=devnull,stderr=devnull))
            end
        end
    end
end

@testset "raw ROOT TH1 options, overlays, and attached functions" begin
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            h=Hist1D([0.,1.,2.,4.],[2.,0.,4.],[2.,0.,5.])
            spec=GraphSpec(legend=:none)
            files=Dict{String,String}()
            for option in ("", "PIE", "HBAR3", "BAR2", "AXIS", "HIST ][",
                           "E0 E1 X0", "P0", "LF2", "TEXT45")
                key=isempty(option) ? "default" : replace(option,r"[^A-Za-z0-9]"=>"_")
                output=render_root(h,spec,joinpath(dir,key);draw=option)
                @test isfile(output.png) && filesize(output.png)>1000
                files[option]=output.png
            end
            @test read(files["PIE"]) != read(files["HBAR3"])
            @test read(files["AXIS"]) != read(files["P0"])
            loghist=Hist1D([1.,2.,3.],[0.,5.],[0.,5.])
            @test isfile(render_root(loghist,GraphSpec(logx=true,logy=true,legend=:none),
                                      joinpath(dir,"raw-log");draw="E1").png)
            other=Hist1D([0.,1.,2.,4.],[1.,3.,2.],[1.,3.,2.])
            over=render_root([h,other],spec,joinpath(dir,"overlay");
                             draw=["HIST","E1 SAME"],
                             styles=[GraphStyle(color="#336699"),GraphStyle(color="#cc3311")])
            @test isfile(over.png) && filesize(over.png)>1000
            @test_throws ArgumentError render_root([h,other],spec,joinpath(dir,"bad-overlay");
                                                   draw=["HIST","E1"])
            fit=fit_root(Hist1D([0.,1.,2.,3.],[5.,6.,4.],[5.,6.,4.]),:pol0;
                         initial=[5.])
            @test isfile(render_root(h,spec,joinpath(dir,"func");
                                      draw="FUNC",attached_fit=fit).png)
            exe=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            input=joinpath(dir,"raw.tsv")
            ROOTPlots.write_hist_draw_spec(input,h,spec,
                ROOTPlots.raw_hist_draw_limits(h,spec,GraphStyle(),"H"),
                "PIE","H",GraphStyle())
            valid=read(input,String)
            @test occursin("DRAWHISTRAW\t",valid)
            for malformed in (replace(valid,"\tPIE\t"=>"\t"*repeat("A",128)*"\t"),
                              replace(valid,"\tPIE\t"=>"\tPIE SAME\t"),
                              replace(valid,"\tPIE\t"=>"\tSPEC BF(1600\t"),
                              replace(valid,"\tPIE\t"=>"\tSPEC \\\";evil()\t"))
                write(input,malformed)
                @test !success(pipeline(Cmd([exe,input,joinpath(dir,"invalid")]);
                                        stdout=devnull,stderr=devnull))
            end
        end
    end
end

@testset "native ROOT TH1 draw modes" begin
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            h = Hist1D([0.,1,2,4], [2.,3,1], [2.,3,1])
            spec = GraphSpec(xlabel="x",ylabel="Count",legend=:none)
            outputs = Dict{Symbol,String}()
            for mode in (:hist,:e0,:e1,:e2,:e3,:e4,:e5,:e6,:hist_e1,
                         :bar,:points,:points_zero,:line,:curve,:text)
                files = render_root(h,spec,joinpath(dir,String(mode));draw=mode)
                @test isfile(files.png) && filesize(files.png) > 1000
                outputs[mode] = files.png
            end
            @test read(outputs[:hist]) != read(outputs[:e1])
            @test read(outputs[:e1]) != read(outputs[:e2])
            @test read(outputs[:hist]) != read(outputs[:bar])
            pdf = render_root(h,spec,joinpath(dir,"pdf");draw="HIST E1",format=:pdf,
                              style=GraphStyle(color="#cc3311"))
            @test pdf.png === nothing && isfile(pdf.pdf)
            @test read(pdf.pdf,4) == collect(codeunits("%PDF"))
            z = Hist1D([1.,2,3,4], [0.,5,2], [0.,5,2])
            logfiles = render_root(z,GraphSpec(logx=true,logy=true,legend=:none),
                                   joinpath(dir,"log-zero");draw=:e1)
            @test isfile(logfiles.png) && filesize(logfiles.png) > 1000
            legacy = render_root(h,spec,joinpath(dir,"steps"))
            explicit = render_root(h,spec,joinpath(dir,"explicit-steps");draw=:steps)
            @test read(legacy.png) == read(explicit.png)
        end
    end
end

@testset "native graph protocol rejects malformed records" begin
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            exe = ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            src = joinpath(dir,"valid.tsv")
            s = GraphSeries([1.0,2.0],[2.0,3.0])
            spec = GraphSpec()
            ROOTPlots.write_graph_spec(src,[s],spec,ROOTPlots.graph_limits([s],spec))
            valid = read(src,String)
            for malformed in (
                replace(valid,"DRAW\t800\t600\t0\t0"=>"DRAW\t800\t600\t2\t0"),
                replace(valid,"POINT\t1\t1"=>"POINT\t1\t2"),
                replace(valid,"END\t1"=>"UNKNOWN\t1"),
                valid * "END\t1\n",
            )
                path = joinpath(dir,"bad.tsv")
                write(path,malformed)
                cmd = Cmd([exe,path,joinpath(dir,"bad-output")])
                @test !success(pipeline(cmd;stdout=devnull,stderr=devnull))
            end
        end
    end
end

@testset "ROOT fit result booleans are strict" begin
    mktempdir() do dir
        path=joinpath(dir,"result.tsv")
        valid="ROOTPlotsFitResult\t1\nFIT\t0\t1\t3\t1.0\t2\t3\t0.0\tMinuit2\t0\n" *
              "PARAM\t1\t2.0\t0.1\t0\nCOV\t1\t0.01\nEND\n"
        write(path,valid)
        @test ROOTPlots.read_fit_result(path,:pol0,(0.,1.),:chi2).valid
        for bad in (replace(valid,"FIT\t0\t1"=>"FIT\t0\tgarbage"),
                    replace(valid,"0.1\t0"=>"0.1\tgarbage"),
                    replace(valid,"Minuit2\t0"=>"Minuit2\tgarbage"))
            write(path,bad)
            @test_throws Exception ROOTPlots.read_fit_result(path,:pol0,(0.,1.),:chi2)
        end
    end
end

@testset "2-D grids and 3-D XYZ inputs" begin
    xedges=[0.,1.,3.]; yedges=[-2.,0.,4.,7.]
    values=[1. 2. 3.; 4. 5. 6.]
    grid=Grid2D(xedges,yedges,values)
    xedges[1]=-100; values[1,1]=-100
    @test grid.xedges[1]==0 && grid.values[1,1]==1
    @test size(grid.values)==(2,3)
    @test_throws ArgumentError Grid2D([0.,1.],[0.,1.],values)
    @test_throws ArgumentError Grid2D([0.,1.,1.],[0.,1.,2.,3.],values)
    @test_throws ArgumentError Grid2D([0.,1.,3.],[0.,1.,2.,3.],fill(NaN,2,3))
    @test_throws ArgumentError Graph3D([0.,1.],[0.],[0.,1.])
    @test_throws ArgumentError Graph3D([0.],[0.],[Inf])
    @test_throws ArgumentError Graph3D([1.],[2.],[3.])
    @test_throws ArgumentError Graph3D([1.,1.],[0.,2.],[0.,1.])
    @test_throws ArgumentError ROOTPlots.validate_nd_spec(Graph2DSpec(draw=:bad))
    @test_throws ArgumentError ROOTPlots.validate_nd_spec(Graph2DSpec(title="bad\tlabel"))
    @test_throws ArgumentError ROOTPlots.validate_nd_spec(Graph3DSpec(theta=Inf))
    mktemp() do path, io
        close(io)
        ROOTPlots.write_grid_spec(path,grid,Graph2DSpec())
        spec=read(path,String)
        @test occursin("GRID\t2\t3\n",spec)
        @test occursin("CELL\t1\t2\t2.0\n",spec)
        @test occursin("CELL\t2\t1\t4.0\n",spec)
    end
    grid.values[1,1]=NaN
    @test_throws ArgumentError ROOTPlots.validate_grid(grid)
end

@testset "3-D binned histogram ownership and validation" begin
    edges=[0.,1.,2.]
    values=reshape(collect(1.0:8.0),2,2,2)
    grid=Grid3D(edges,edges,edges,values)
    edges[1]=-1; values[1,1,1]=-1
    @test grid.xedges[1]==0 && grid.values[1,1,1]==1
    @test_throws ArgumentError Grid3D([0.,1.],[0.,1.,2.],[0.,1.,2.],ones(2,2,2))
    @test_throws ArgumentError Grid3D([0.,1.,1.],[0.,1.,2.],[0.,1.,2.],ones(2,2,2))
    @test_throws ArgumentError Grid3D([0.,1.,2.],[0.,1.,2.],[0.,1.,2.],fill(NaN,2,2,2))
    @test_throws ArgumentError render_root(grid,Graph3DSpec(),
                                           joinpath(tempdir(),"invalid-volume"))
end

@testset "raw TH2 and TH3 option families" begin
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            e=[0.,1.,2.,3.,4.]
            a=Grid2D(e,e,[exp(-((i-2)^2+(j-2)^2)) for i=1:4,j=1:4])
            b=Grid2D(e,e,[exp(-((i-3)^2+(j-3)^2)) for i=1:4,j=1:4])
            for option in ("COLZ2", "BOX1", "CONT4", "SURF2Z", "LEGO2",
                           "ARR", "TEXT45", "AITOFF", "CANDLEX", "VIOLINY")
                output=render_root(a,Graph2DSpec(draw=option),
                                   joinpath(dir,"th2-"*lowercase(option)))
                @test isfile(output.png) && filesize(output.png)>1000
            end
            overlay=render_root([a,b],Graph2DSpec(contours=5),joinpath(dir,"th2-overlay");
                                draw=["COLZ","CONT3 SAME"],
                                styles=[GraphStyle(),GraphStyle(color="#cc3311")])
            @test isfile(overlay.png) && filesize(overlay.png)>1000
            @test_throws ArgumentError render_root(a,Graph2DSpec(draw="COLZ SAME"),
                                                   joinpath(dir,"bad-grid"))
            @test_throws ArgumentError render_root(a,Graph2DSpec(draw="GLLEGO"),
                                                   joinpath(dir,"bad-gl-grid"))
            volume=Grid3D(e,e,e,[exp(-0.5*((i-2.5)^2+(j-2.5)^2+(k-2.5)^2))
                                    for i=1:4,j=1:4,k=1:4])
            for option in ("BOX", "BOX1", "BOX2Z", "BOX3", "ISO", "LEGO")
                output=render_root(volume,Graph3DSpec(draw=option),
                                   joinpath(dir,"th3-"*lowercase(option)))
                @test isfile(output.png) && filesize(output.png)>1000
            end
            @test_throws ArgumentError render_root(volume,Graph3DSpec(draw="ISO SAME"),
                                                   joinpath(dir,"bad-volume"))
            exe=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            input=joinpath(dir,"volume.tsv")
            ROOTPlots.write_volume_spec(input,volume,Graph3DSpec(draw="BOX2Z"))
            valid=read(input,String)
            @test occursin("VOLUME\t4\t4\t4",valid)
            for malformed in (replace(valid,"CELL3\t1\t1\t1"=>"CELL3\t2\t1\t1"),
                              replace(valid,"ZEDGE\t2\t1.0"=>"ZEDGE\t2\t0.0"))
                write(input,malformed)
                @test !success(pipeline(Cmd([exe,input,joinpath(dir,"bad-volume-protocol")]);
                                        stdout=devnull,stderr=devnull))
            end
        end
    end
end

@testset "native ROOT 2-D and 3-D plots" begin
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            grid=Grid2D([0.,1.,3.],[-2.,0.,4.,7.],[0. 2. 3.; 4. -1. 6.])
            cases=((grid,Graph2DSpec(draw=:heatmap)),
                   (grid,Graph2DSpec(draw=:contour,contours=5)),
                   (grid,Graph3DSpec(draw=:surface)),
                   (Graph3D([0.,1.,1.,2.],[0.,1.,1.,2.],[0.,1.,2.,3.]),
                    Graph3DSpec(draw=:points)))
            for (i,(data,spec)) in enumerate(cases)
                output=render_root(data,spec,joinpath(dir,"plot$i"))
                @test isfile(output.png) && filesize(output.png)>1000
                @test output.pdf === nothing
                @test !isfile(joinpath(dir,"plot$i.pdf"))
            end
            @test isfile(render_root(Grid2D([0.,1.],[0.,1.],zeros(1,1)),
                                          Graph2DSpec(),joinpath(dir,"zero")).png)
            @test_throws ArgumentError render_root(Grid2D([0.,1.],[0.,1.],ones(1,1)),
                                                   Graph3DSpec(),joinpath(dir,"one-cell-surface"))
            @test_throws ArgumentError render_root(grid,Graph3DSpec(draw=:points),joinpath(dir,"bad"))
            @test_throws ArgumentError render_root(last(cases)[1],Graph3DSpec(),joinpath(dir,"bad"))

            exe=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            path=joinpath(dir,"grid.tsv")
            ROOTPlots.write_grid_spec(path,grid,Graph2DSpec())
            valid=read(path,String)
            for malformed in (replace(valid,"CELL\t1\t2"=>"CELL\t2\t2"),
                              replace(valid,"GRID\t2\t3"=>"GRID\t2\t4"),
                              valid * "END\t1\n")
                write(path,malformed)
                @test !success(pipeline(Cmd([exe,path,joinpath(dir,"invalid")]);
                                        stdout=devnull,stderr=devnull))
            end
        end
    end
end

@testset "select PNG or PDF output" begin
    mktempdir() do dir
        grid = Grid2D([0.,1.,2.],[0.,1.,2.],[1. 2.; 3. 4.])
        cloud = Graph3D([0.,1.,0.5],[0.,0.5,1.],[1.,2.,3.])
        graph = GraphSeries([0.,1.],[1.,2.])
        hist = Hist1D([0.,1.,2.],[2.,3.],[2.,3.])
        cmp = comparison(hist,[Component("MC",hist,"#123456")])
        badstem = joinpath(dir,"not-created","plot")
        for request in (() -> render_root(grid,Graph2DSpec(),badstem;
                                          format=:svg,root_config="/nonexistent"),
                        () -> render_root([graph],GraphSpec(),badstem;
                                          format=:svg,root_config="/nonexistent"),
                        () -> render_root([cmp],[PanelSpec("x")],badstem;
                                          format=:svg,root_config="/nonexistent"))
            @test_throws ArgumentError request()
        end
        @test !ispath(dirname(badstem))

        if TEST_NATIVE_ROOT
            cases = (([cmp],[PanelSpec("x")]),
                     ([graph],GraphSpec()),
                     (grid,Graph2DSpec(draw=:heatmap)),
                     (grid,Graph2DSpec(draw=:contour)),
                     (grid,Graph3DSpec(draw=:surface)),
                     (cloud,Graph3DSpec(draw=:points)))
            default = render_root([cmp],[PanelSpec("x")],joinpath(dir,"root-default"))
            @test isfile(default.png) && default.pdf === nothing
            @test !isfile(joinpath(dir,"root-default.pdf"))
            both = render_root([cmp],[PanelSpec("x")],joinpath(dir,"root-both");format=:both)
            @test isfile(both.png) && isfile(both.pdf)
            for (i,(data,spec)) in enumerate(cases), format in (:png,:pdf)
                stem = joinpath(dir,"root-$(i)-$(format)")
                result = render_root(data,spec,stem;format)
                selected = getproperty(result,format)
                other = format == :png ? :pdf : :png
                @test selected == stem * ".$(format)"
                @test isfile(selected) && filesize(selected) > 1000
                @test getproperty(result,other) === nothing
                @test !isfile(stem * ".$(other)")
                magic = format == :png ? UInt8[0x89,0x50,0x4e,0x47] : collect(codeunits("%PDF-"))
                @test read(selected)[1:length(magic)] == magic
            end
            stem = joinpath(dir,"preserve-root")
            write(stem * ".pdf","older PDF")
            result = render_root(grid,Graph2DSpec(),stem;format=:png)
            @test result.pdf === nothing && read(stem * ".pdf",String) == "older PDF"
        end

    end
end

@testset "native ROOT graph drawing and fitting" begin
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            x = collect(0.0:4.0)
            y = [1.1,2.9,5.2,6.8,9.1]
            g = GraphSeries(x,y;yerrors=fill(0.2,5),label="Data")
            f = fit_root(g,:pol1;initial=[0.0,1.0])
            @test f.status == 0 && f.valid && f.covariance_status == 3
            @test f.parameters ≈ [1.04,1.99] atol=1e-5
            @test f.ndf == 3 && f.npoints == 5
            @test size(f.covariance) == (2,2) && f.covariance[1,1] > 0
            @test !isempty(f.minimizer) # implementation name can vary across ROOT releases
            bounded = fit_root(g,:pol1;initial=[1.0,0.2],
                               bounds=[nothing,(0.0,0.5)])
            @test bounded.status == 0 && bounded.parameters[2] <= 0.500001

            gx=collect(-3.0:0.5:3.0)
            gy=10 .* exp.(-0.5 .* gx.^2)
            gaussian_graph=GraphSeries(gx,gy;yerrors=fill(0.2,length(gx)))
            gaussian_hist=Hist1D(collect(-3.25:0.5:3.25),gy,fill(0.04,length(gy)))
            for data in (gaussian_graph,gaussian_hist)
                fixed_gauss=fit_root(data,:gaus;initial=[10.,0.,2.],fixed=[3])
                @test fixed_gauss.status==0 && fixed_gauss.parameters[3]==2 &&
                      fixed_gauss.fixed[3] && all(iszero,fixed_gauss.covariance[3,:])
                bounded_gauss=fit_root(data,:gaus;initial=[10.,0.,2.],
                                       bounds=[nothing,nothing,(1.5,2.5)])
                @test bounded_gauss.status==0 && 1.499999 <= bounded_gauss.parameters[3] <= 2.5
            end

            files = render_root([fit_curve(f),g],GraphSpec(xlabel="x",ylabel="y",
                                 title="Graph fit"),joinpath(dir,"fit"))
            @test isfile(files.png) && filesize(files.png)>1000
            @test files.pdf === nothing && !isfile(joinpath(dir,"fit.pdf"))
            asym = GraphSeries([1.0,2,3],[2.0,4,8];
                               xerrors=(low=[0.1,0.1,0.1],high=[0.2,0.2,0.2]),
                               yerrors=(low=[0.2,0.4,0.8],high=[0.3,0.6,1.2]),
                               label="Asymmetric")
            af = fit_root(asym,:pol1;initial=[0.0,2.0],xerrors=:ignore)
            @test af.status == 0 && af.npoints == 3
            logfiles = render_root([asym],GraphSpec(logx=true,logy=true,
                                   legend=:topleft),joinpath(dir,"log"))
            @test isfile(logfiles.png) && logfiles.pdf === nothing
            @test !isfile(joinpath(dir,"log.pdf"))
            visible=GraphSeries([1.,2.],[2.,4.];yerrors=[0.5,0.5],errors=:bars)
            hidden=GraphSeries([1.,2.],[2.,4.];yerrors=[0.5,0.5],errors=:none)
            axes=GraphSpec(xlimits=(0.5,2.5),ylimits=(1.,5.),legend=:none)
            shown=render_root([visible],axes,joinpath(dir,"bars"))
            suppressed=render_root([hidden],axes,joinpath(dir,"no-bars"))
            @test read(shown.png) != read(suppressed.png)
            crossing=GraphSeries([1.,2.],[1.,2.];yerrors=[2.,2.],errors=:none)
            @test isfile(render_root([crossing],GraphSpec(logy=true,legend=:none),
                                      joinpath(dir,"log-hidden")).png)

            hist = Hist1D([0.,1,2,3,4],[2.,0,4,2],[2.,0,4,2])
            hf = fit_root(hist,:pol0;initial=[1.0],objective=:poisson)
            @test hf.status == 0 && hf.npoints == 4 && hf.ndf == 3
            @test hf.parameters[1] ≈ 2.0 atol=1e-3
            @test hf.chi2_kind == :poisson_deviance
            hc = fit_root(hist,:pol0;initial=[1.0])
            @test hc.npoints == 3 # ROOT chi2 excludes the zero-error bin
            pulls = fit_residuals(hist,hc)
            @test count(pulls.valid) == hc.npoints
            @test sum(abs2, pulls.pull[pulls.valid]) ≈ hc.chi2 atol=1e-6
            hist_graph = histogram_points(hist)
            hist_line = histogram_steps(hist)
            pull_graph = fit_pull_points(hist, hc)
            @test length(pull_graph.x) == hc.npoints
            @test sum(abs2, pull_graph.y) ≈ hc.chi2 atol=1e-6
            overlay = render_root([fit_curve(hc), hist_graph],
                                  GraphSpec(xlabel="x", ylabel="Bin count"),
                                  joinpath(dir,"hist-fit"))
            pull_plot = render_root([pull_graph],
                                    GraphSpec(xlabel="x", ylabel="Pull"),
                                    joinpath(dir,"hist-pull"))
            @test isfile(overlay.png) && filesize(overlay.png) > 1000
            standalone = render_root(hist, GraphSpec(xlabel="x", ylabel="Bin count"),
                                     joinpath(dir,"hist-steps"))
            @test isfile(standalone.png) && filesize(standalone.png) > 1000
            stepped_overlay = render_root([hist_line, hist_graph],
                                          GraphSpec(xlabel="x", ylabel="Bin count"),
                                          joinpath(dir,"hist-steps-errors"))
            @test isfile(stepped_overlay.png) && filesize(stepped_overlay.png) > 1000
            @test isfile(pull_plot.png) && filesize(pull_plot.png) > 1000
            ranged = fit_root(hist,:pol0;initial=[1.0],range=(0.2,2.6))
            ranged_pulls = fit_residuals(hist,ranged)
            @test ranged_pulls.valid == [true,false,true,false]
            @test count(ranged_pulls.valid) == ranged.npoints
            @test sum(abs2, ranged_pulls.pull[ranged_pulls.valid]) ≈ ranged.chi2 atol=1e-6
            combined = render_root(hist,ranged,
                                   GraphSpec(xlabel="x",ylabel="Bin count",title="Fit / pull",
                                             xlimits=(0.,4.),legend=:topleft),
                                   joinpath(dir,"fit-pull");format=:both)
            @test isfile(combined.png) && filesize(combined.png) > 1000
            @test isfile(combined.pdf) && read(combined.pdf,4) == collect(codeunits("%PDF"))
            @test fit_pull_points(hist,ranged).x == [0.5,2.5]
            pdf_only = render_root(hist,ranged,GraphSpec(),joinpath(dir,"fit-pull-pdf");
                                   format=:pdf)
            @test pdf_only.png === nothing && isfile(pdf_only.pdf)
            @test !isfile(joinpath(dir,"fit-pull-pdf.png"))
            @test_throws ArgumentError render_root(hist,hf,GraphSpec(),joinpath(dir,"bad-poisson"))
            positive = Hist1D([1.,2,4,8],[4.,6,5],[4.,6,5])
            positive_fit = fit_root(positive,:pol0;initial=[5.])
            log_combined = render_root(positive,positive_fit,
                                       GraphSpec(logx=true,logy=true,legend=:none),
                                       joinpath(dir,"fit-pull-log"))
            @test isfile(log_combined.png) && filesize(log_combined.png) > 1000
            @test_throws ArgumentError render_root(hist,ranged,GraphSpec(logx=true),
                                                   joinpath(dir,"bad-logx"))
            @test 0 <= chi2_pvalue(hc) <= 1
            @test 0 <= chi2_pvalue(hf) <= 1
            poisson_deviance = 2sum(n == 0 ? hf.parameters[1] :
                                    hf.parameters[1]-n+n*log(n/hf.parameters[1])
                                    for n in hist.sumw)
            @test hf.chi2 ≈ poisson_deviance atol=1e-8
            weighted = Hist1D([0.,1,2,3],[1.5,2.5,3.5],[1.1,1.8,2.9])
            wf = fit_root(weighted,:pol0;initial=[2.0],objective=:weighted_poisson)
            @test wf.status == 0 && wf.npoints == 3
            @test wf.chi2_kind == :weighted_poisson_objective
            @test_throws ArgumentError chi2_pvalue(wf)
            plain = fit_root(GraphSeries(x,y),:pol1;initial=[0.,1.])
            @test plain.normalized_errors
            @test_throws ArgumentError chi2_pvalue(plain)
            @test_throws ArgumentError fit_root(Hist1D([0.,1,2],[2.,0],[0.,0]),
                                                :pol0;initial=[1.],objective=:weighted_poisson)
            @test_throws ArgumentError fit_root(weighted,:pol0;initial=[2.0],objective=:poisson)
            @test_throws ArgumentError fit_root(g,:bad;initial=[1])
            @test_throws ArgumentError fit_root(GraphSeries(x,y;yerrors=zeros(5)),
                                               :pol1;initial=[0.,1.])
            @test_throws ArgumentError fit_root(hist,:pol0;initial=[1.0],
                                                 range=(0.0,0.0))
        end
    end
end

@testset "threaded DataFrame filling is reproducible and validates events" begin
    n = 200_000
    df = DataFrame(x=[mod(i, 100) / 10 for i in 1:n],
                   w=[mod(i, 7) == 0 ? -0.25 : 0.5 for i in 1:n])
    edges = [0.0, 1.0, 3.0, 10.0]
    off = histogram(df, :x, edges; weight=:w, threads=:off)
    on = histogram(df, :x, edges; weight=:w, threads=:on)
    auto = histogram(df, :x, edges; weight=:w)
    @test isequal(off.hist.sumw, on.hist.sumw) && isequal(on.hist.sumw, auto.hist.sumw)
    @test isequal(off.hist.sumw2, on.hist.sumw2) && isequal(on.hist.sumw2, auto.hist.sumw2)
    @test off.underflow == on.underflow == auto.underflow
    @test off.overflow == on.overflow == auto.overflow
    subset_on = histogram(view(df, 1:100, :), :x, edges; weight=:w, threads=:on)
    subset_off = histogram(df[1:100, :], :x, edges; weight=:w, threads=:off)
    @test subset_on.hist.sumw == subset_off.hist.sumw
    @test histogram(DataFrame(x=Float64[]), :x, edges; threads=:on).hist.sumw == zeros(3)
    @test_throws ArgumentError histogram(df, :x, edges; threads=:bad)
    @test_throws ArgumentError histogram(df, :missing, edges)

    bad = DataFrame(x=Union{Missing,Float64}[1.0, missing, 2.0], w=[1.0, 1.0, 1.0])
    @test occursin("row 2", sprint(showerror,
        try histogram(bad, :x, edges; weight=:w, threads=:on)
        catch err; err end))
    multi_bad = DataFrame(x=Vector{Union{Missing,Float64}}(df.x))
    multi_bad.x[70_000] = missing
    multi_bad.x[150_000] = missing
    failure = try histogram(multi_bad, :x, edges; threads=:on)
    catch err; err end
    @test failure isa ArgumentError && occursin("row 70000", sprint(showerror, failure))
    @test_throws ArgumentError histogram(DataFrame(x=[1.0], w=[1e200]), :x, edges; weight=:w)
    huge = histogram(DataFrame(x=[1.0], w=[typemax(Int)]), :x, edges; weight=:w)
    @test huge.hist.sumw2[2] == Float64(typemax(Int))^2
    @test_throws ArgumentError histogram(DataFrame(x=[big"1e400"]), :x, edges)
end

@testset "stack sums, replacement uncertainty, and ratio" begin
    edges = [0.0, 1.0, 2.0, 4.0]
    data = Hist1D(edges, [12, 0, 4], [12, 0, 4])
    a = Component("A", Hist1D(edges, [7, 0, 1], [7, 0, 1]), "#123456")
    b = Component("B", Hist1D(edges, [3, 0, 3], [3, 0, 3]), "#654321")
    cmp = comparison(data, [a, b]; data_errors=BinErrors([3, 0, 1], [4, 1.84, 2]),
                     prediction_errors=BinErrors([1, 0, 2], [2, 0, 3]))
    @test cmp.cumulative == [[7, 0, 1], [10, 0, 4]]
    @test cmp.prediction == [10, 0, 4]
    @test cmp.prediction_errors.low == [1, 0, 2] # no double counting
    @test cmp.valid == [true, false, true]
    @test cmp.ratio[1] ≈ 1.2
    @test isnan(cmp.ratio[2])
    @test cmp.ratio[3] ≈ 1.0
    @test cmp.ratio_data_errors.low[1] ≈ 0.3
    @test cmp.ratio_data_errors.high[1] ≈ 0.4
    @test cmp.ratio_prediction_errors.high[3] ≈ 0.75
    @test_throws ArgumentError comparison(data, [Component("bad", Hist1D([0, 1], [1], [1]), "#abcdef")])
    @test_throws ArgumentError Component("signed", Hist1D([0, 1], [-1], [1]), "#abcdef")

    zero_data = Hist1D([0.0, 1.0], [0.0], [0.0])
    positive = Component("MC", Hist1D([0.0, 1.0], [2.0], [1.0]), "#abcdef")
    zero_cmp = comparison(zero_data, [positive]; data_errors=BinErrors([0.0], [1.84]))
    @test zero_cmp.valid == [true]
    @test zero_cmp.ratio == [0.0]
    @test zero_cmp.ratio_data_errors.high == [0.92]
end

@testset "native ROOT is optional; no remote fallback" begin
    @test root_available() isa Bool
    if !root_available()
        data = Hist1D([0, 1], [1], [1])
        component = Component("MC", Hist1D([0, 1], [1], [1]), "#abcdef")
        @test_throws ArgumentError render_root([comparison(data, [component])],
                                              [PanelSpec("x")], tempname())
    end
end

@testset "ROOT renderer specification preserves labels and numeric settings" begin
    data = Hist1D([0.0, 1.0], [2.0], [2.0])
    component = Component("Model", Hist1D([0.0, 1.0], [2.0], [1.0]), "#abcdef")
    cmp = comparison(data, [component]; uncertainty_label="Illustrative unc.")
    mktemp() do path, io
        close(io)
        ROOTPlots.write_root_spec(path, [cmp], [PanelSpec("m(eμ)"; annotation="eμ channel")],
                                  (1//2, 3//2))
        spec = read(path, String)
        @test occursin("SETTINGS\t0.5\t1.5", spec)
        @test startswith(spec, "ROOTPlots\t2\n")
        @test occursin("STYLE\t20\t0.8\t#000000\t1\t1001\t3345\t#000000\t2", spec)
        @test occursin("AXIS\t1\t0\t0\t0", spec)
        @test occursin("PANEL\t1\tm(e#mu)\tEvents\te#mu channel\t\tIllustrative unc.\t1", spec)
    end
end

@testset "ROOT-style panel options and log limits" begin
    edges = [1.0, 2.0, 4.0]
    data = Hist1D(edges, [0.0, 3.0], [0.0, 3.0])
    mc = Component("MC", Hist1D(edges, [0.0, 2.0], [0.0, 2.0]), "#abcdef")
    cmp = comparison(data, [mc]; data_errors=BinErrors([0.0, 1.0], [1.84, 1.0]),
                     prediction_errors=BinErrors([0.0, 3.0], [1.0, 3.0]))
    spec = PanelSpec("x"; logx=true, logy=true, ratio_logy=true,
                     ymin=0.1, ratio_limits=(0.1, 2.0))
    lim = ROOTPlots.panel_limits(cmp, spec, (0.5, 1.5))
    @test lim.ymin == 0.1 && lim.ratio_min == 0.1 && lim.ratio_max == 2.0
    @test_throws ArgumentError ROOTPlots.panel_limits(cmp, PanelSpec("x"; logy=true, ymin=0), (0.5,1.5))
    @test_throws ArgumentError ROOTPlots.panel_limits(cmp, PanelSpec("x"; ratio_logy=true, ratio_limits=(0,2)), (0.5,1.5))
    zero_edge = Hist1D([0.0, 2.0], [1.0], [1.0])
    zero_cmp = comparison(zero_edge, [Component("MC", zero_edge, "#abcdef")])
    @test_throws ArgumentError ROOTPlots.panel_limits(zero_cmp, PanelSpec("x"; logx=true), (0.5,1.5))
    @test ROOTPlots.validate_style(PlotStyle(marker_style=19)).marker_style == 19
    mktemp() do path, io
        close(io)
        ROOTPlots.write_root_spec(path, [cmp], [spec], (0.5,1.5),
                                  PlotStyle(marker_style=21, line_width=2))
        text = read(path, String)
        @test occursin("STYLE\t21\t0.8\t#000000\t2", text)
        @test occursin("AXIS\t1\t1\t1\t1\t0.1", text)
    end
end

end # ROOTPlots
