@testset "graph presets, fills, and callable sampling" begin
    legacy=GraphStyle("#123456",21,0.9,3,2)
    @test legacy.color=="#123456" && legacy.fill_color===nothing
    @test ROOTPlots.check_graph_style(legacy)===legacy
    for invalid in (GraphStyle(fill_color="red"),GraphStyle(fill_alpha=NaN),
                    GraphStyle(fill_alpha=1.1),GraphStyle(fill_style=-1))
        @test_throws ArgumentError ROOTPlots.check_graph_style(invalid)
    end
    x=[1.,2.,3.,4.]; y=[2.,4.,3.,5.]
    for draw in (:points,:line,:linepoints,:curve,:curvepoints,:fill,:bar,:bar_offset,"C","P3")
        @test ROOTPlots.validate_graph(GraphSeries(x,y;draw)) isa GraphSeries
    end
    for errors in (:bars,:none,:band,:smoothband,:boxes,:caps,:nocaps,:arrows,:filledarrows)
        s=GraphSeries(x,y;yerrors=(low=fill(.2,4),high=fill(.4,4)),errors)
        @test ROOTPlots.graph_limits([s],GraphSpec())[2][2]>maximum(y)
    end
    @test_throws ArgumentError GraphSeries(x,y;errors=:band)
    for option in ("AP","P SAME","P; system()","P\n","P\0")
        @test_throws ArgumentError GraphSeries(x,y;draw=option)
    end
    low=GraphSeries([1.,2.],[.2,.3];yerrors=[.3,.1],errors=:band)
    @test_throws ArgumentError ROOTPlots.graph_limits([low],GraphSpec(logy=true))
    @test ROOTPlots.graph_limits([GraphSeries(x,y;draw=:bar)],GraphSpec())[2][1]<0
    sampled=function_curve(t->t^2,(1.,4.);n=7,label="Julia sampled")
    @test sampled.x==collect(range(1.,4.;length=7))
    @test sampled.y==sampled.x.^2 && sampled.draw==:line && sampled.errors==:none
    @test_throws ArgumentError function_curve(identity,(0.,1.);n=1)
    @test_throws ArgumentError function_curve(identity,(1.,0.))
    @test_throws ArgumentError function_curve(t->Inf,(0.,1.))
    @test_throws ArgumentError function_curve(t->complex(t,1),(0.,1.))
    p=[1.,2.]; model=ModelCurve(:pol1,p;range=(1.,4.))
    p[1]=100.
    @test model.parameters==[1.,2.]
    @test_throws ArgumentError ModelCurve(:unknown,[1.];range=(0.,1.))
    @test_throws ArgumentError ModelCurve(:gaus,[1.,0.,0.];range=(-1.,1.))
    @test_throws ArgumentError ModelCurve(:pol1,[1.,2.];range=(0.,1.),npx=1)
    mktempdir() do dir
        plain=GraphSeries(x,y)
        path=joinpath(dir,"legacy.tsv")
        ROOTPlots.write_graph_spec(path,[plain],GraphSpec(),ROOTPlots.graph_limits([plain],GraphSpec()))
        row=only(filter(l->startswith(l,"SERIES\t"),readlines(path)))
        @test length(split(row,'\t'))==12
        styled=GraphSeries(x,y;yerrors=fill(.3,4),errors=:band,
                           style=GraphStyle(fill_color="#336699",fill_alpha=.6,fill_style=3004))
        ROOTPlots.write_graph_spec(path,[styled],GraphSpec(),ROOTPlots.graph_limits([styled],GraphSpec()))
        row=only(filter(l->startswith(l,"SERIES\t"),readlines(path)))
        @test length(split(row,'\t'))==15 && endswith(row,"#336699\t3004\t0.6")
        if TEST_NATIVE_ROOT
            spec=GraphSpec(xlabel="x",ylabel="y",legend=:none)
            for (i,errors) in enumerate((:bars,:band,:smoothband,:boxes,:caps,:nocaps,:arrows,:filledarrows))
                g=GraphSeries(x,y;yerrors=fill(.3,4),errors,style=GraphStyle(color="#336699"))
                files=render_root([g],spec,joinpath(dir,"errors-$i");format=:both)
                @test read(files.png,8)==UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a]
                @test read(files.pdf,4)==collect(codeunits("%PDF"))
            end
            for (i,draw) in enumerate((:curve,:curvepoints,:fill,:bar,:bar_offset,"C"))
                files=render_root([GraphSeries(x,y;draw)],spec,joinpath(dir,"draw-$i"))
                @test filesize(files.png)>1000
            end
            files=render_root(model,spec,joinpath(dir,"native-tf1");format=:both)
            @test filesize(files.png)>1000 && read(files.pdf,4)==collect(codeunits("%PDF"))
            executable=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            ROOTPlots.write_graph_spec(path,[styled],spec,ROOTPlots.graph_limits([styled],spec))
            original=read(path,String)
            for broken in (replace(original,"#336699\t3004\t0.6"=>"#336699\t3004\t2.0"),
                           replace(original,"points\tband"=>"ROOT:AP\tband"))
                write(path,broken)
                @test !success(pipeline(Cmd([executable,path,joinpath(dir,"bad")]);stdout=devnull,stderr=devnull))
            end
        end
    end
end

@testset "scattered native XYZ interpolation and errors" begin
    x=[0.,1.,0.,1.,.4]; y=[0.,0.,1.,1.,.6]; z=1 .+ x.^2 .+ y
    cloud=Graph3D(x,y,z)
    errors=Graph3DErrors(x,y,z;xerrors=fill(.02,5),yerrors=fill(.03,5),zerrors=fill(.1,5))
    x[1]=-99.
    @test errors.x[1]==0.
    @test_throws ArgumentError Graph3DErrors(cloud.x,cloud.y,cloud.z;zerrors=[.1])
    @test_throws ArgumentError Graph3DErrors(cloud.x,cloud.y,cloud.z;zerrors=fill(-.1,5))
    @test_throws ArgumentError ROOTPlots.validate_xyz_plot(cloud,Graph3DSpec(draw=:surface))
    @test_throws ArgumentError ROOTPlots.validate_xyz_plot(cloud,Graph3DSpec(draw=:points,logz=true))
    @test_throws ArgumentError ROOTPlots.validate_xyz_plot(cloud,Graph3DSpec(draw="ERR"))
    @test_throws ArgumentError ROOTPlots.validate_xyz_plot(cloud,Graph3DSpec(draw="BOX"))
    @test_throws ArgumentError ROOTPlots.validate_nd_spec(Graph3DSpec(resolution=(3,40)))
    duplicate=Graph3D([0.,1.,0.,0.],[0.,0.,1.,0.],[1.,2.,3.,4.])
    @test_throws ArgumentError ROOTPlots.validate_xyz_plot(duplicate,Graph3DSpec(draw=:surface,interpolate=true))
    collinear=Graph3D([0.,1.,2.],[0.,1.,2.],[1.,2.,1.])
    @test_throws ArgumentError ROOTPlots.validate_xyz_plot(collinear,Graph3DSpec(draw=:surface,interpolate=true))
    @test ROOTPlots.validate_xyz_plot(duplicate,Graph3DSpec(draw=:points))=="P0"
    mktempdir() do dir
        path=joinpath(dir,"xyz.tsv")
        ROOTPlots.write_xyz_spec(path,cloud,Graph3DSpec(draw=:points))
        @test length(split(readlines(path)[2],'\t'))==8
        ROOTPlots.write_xyz_spec(path,errors,Graph3DSpec(draw=:points))
        @test length(split(readlines(path)[2],'\t'))==20
        @test occursin("POINT3E\t1\t0.0",read(path,String))
        if TEST_NATIVE_ROOT
            for mode in (:surface,:triangles,:contour,:heatmap)
                files=render_root(cloud,Graph3DSpec(draw=mode,interpolate=true,resolution=(24,28),
                                  title="Native scattered $mode"),joinpath(dir,String(mode));format=:both)
                @test filesize(files.png)>1000 && read(files.pdf,4)==collect(codeunits("%PDF"))
            end
            file=render_root(errors,Graph3DSpec(draw=:points),joinpath(dir,"xyz-errors");format=:both)
            @test filesize(file.png)>1000 && read(file.pdf,4)==collect(codeunits("%PDF"))
            @test filesize(render_root(cloud,Graph3DSpec(draw=:line),joinpath(dir,"xyz-line")).png)>1000
            ROOTPlots.write_xyz_spec(path,cloud,Graph3DSpec(draw=:surface,interpolate=true))
            original=read(path,String)
            executable=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            for broken in (replace(original,"plain\t1\t2\t1"=>"plain\t0\t2\t1"),
                           replace(original,"POINT3\t5\t0.4\t0.6"=>"POINT3\t5\t0.0\t0.0"))
                write(path,broken)
                @test !success(pipeline(Cmd([executable,path,joinpath(dir,"bad")]);stdout=devnull,stderr=devnull))
            end
        end
    end
end
