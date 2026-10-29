using SpecialFunctions: erf, erfc

@testset "batch composition validates references, styles, and coordinates" begin
    data=GraphSeries([0.,1.,2.],[1.,3.,2.];yerrors=[0.2,0.3,0.2],label="Profile")
    axis=AxisStyle(title_font=42,title_size=0.06,title_offset=1.,label_font=42,
                   label_size=0.05,label_offset=0.015,tick_length=0.03,
                   divisions=505,max_digits=3,center_title=false)
    spec=GraphSpec(xaxis=axis,yaxis=axis,
        canvas=CanvasStyle(left_margin=0.14,right_margin=0.04,bottom_margin=0.14,
                           top_margin=0.05,ticks_x=true,ticks_y=true,
                           frame_line_width=2,end_error_size=0.,line_scale_ps=1.),
        lines=[PlotLine(1.,0.,1.,4.;style=2,label="Nominal")],
        texts=[PlotText(0.58,0.84,"#sigma = 180 #pm 10 #mum")],
        legend_spec=LegendSpec(bounds=(0.18,0.72,0.50,0.90),text_size=0.04,
            entries=[LegendEntry(:series,1,"Profile"),LegendEntry(:line,1,"Nominal";option="l")]))
    @test ROOTPlots.validate_scene(spec;nseries=1) === spec
    @test ROOTPlots.with_graph_spec(spec;ylimits=(0.,5.)).texts === spec.texts
    @test GraphSpec("x","y","",nothing,nothing,false,false,:none,(800,600)).legend == :none
    invalid=(GraphSpec(xaxis=AxisStyle(title_size=-1.)),
        GraphSpec(canvas=CanvasStyle(left_margin=0.6,right_margin=0.5)),
        GraphSpec(texts=[PlotText(1.2,0.5,"outside")]),
        GraphSpec(texts=[PlotText(0.5,0.5,"bad\tlabel")]),
        GraphSpec(texts=[PlotText(0.5,0.5,"bad\0label")]),
        GraphSpec(legend_spec=LegendSpec(entries=[LegendEntry(:line,1,"missing")])),
        GraphSpec(legend_spec=LegendSpec(entries=[LegendEntry(:series,2,"missing")])) )
    mktempdir() do dir
        bad=joinpath(dir,"not-created","plot")
        for s in invalid
            @test_throws ArgumentError render_root([data],s,bad;root_config="/not-root")
        end
        @test !isdir(dirname(bad))
        @test_throws ArgumentError ROOTPlots.graph_limits([GraphSeries([1.,2.],[1.,2.])],
            GraphSpec(logx=true,lines=[PlotLine(0.,1.,1.,2.)]))
        @test ROOTPlots.graph_limits([GraphSeries([1.,2.],[1.,2.])],
            GraphSpec(logx=true,lines=[PlotLine(0.,0.,1.,1.;ndc=true)]))[1][1]>0
        # ROOT precision-3 fonts have pixel, not fractional, sizes.
        @test ROOTPlots.validate_scene(GraphSpec(texts=[PlotText(.2,.8,"pixels";font=43,size=18.)]);nseries=1) isa GraphSpec
        if TEST_NATIVE_ROOT
            files=render_root([data],spec,joinpath(dir,"scene");format=:both)
            @test filesize(files.png)>2000 && read(files.pdf,4)==collect(codeunits("%PDF"))
            plain=render_root([data],GraphSpec(),joinpath(dir,"plain"))
            @test read(files.png)!=read(plain.png)
            h=Hist1D([0.,1.,2.,3.],[2.,3.,1.],[2.,3.,1.])
            @test isfile(render_root(h,spec,joinpath(dir,"styled-hist");draw="E1").png)
            @test isfile(render_root(h,spec,joinpath(dir,"styled-steps")).png)
            curve=GraphSeries([0.,1.,2.,3.],[5.,10.,8.,5.];draw=:line,label="Curve")
            mixedspec=ROOTPlots.with_graph_spec(spec;legend_spec=LegendSpec(entries=[
                LegendEntry(:series,1,"Hist"),LegendEntry(:series,2,"Curve";option="l"),
                LegendEntry(:line,1,"Nominal";option="l")]))
            @test isfile(render_root(h,[curve],mixedspec,joinpath(dir,"mixed");draw="E1").png)
            @test isfile(render_root([h,h],mixedspec,joinpath(dir,"hist-overlay");draw=["HIST","E1 SAME"]).png)
            exe=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            input=joinpath(dir,"scene.tsv")
            ROOTPlots.write_graph_spec(input,[data],spec,ROOTPlots.graph_limits([data],spec))
            original=read(input,String)
            for (i,broken) in enumerate((replace(original,"SCENE\t1"=>"SCENE\t9"),
                replace(original,"LENTRY\tline\t1"=>"LENTRY\tline\t99"),
                replace(original,"ENDSCENE"=>"UNKNOWN"),
                replace(original,"TEXT\t1\t1\t0.58"=>"TEXT\t1\t1\t1.58")))
                write(input,broken)
                @test !success(pipeline(`$exe $input $(joinpath(dir,"bad-$i"))`;stdout=devnull,stderr=devnull))
            end
            fit=fit_root(h,:pol0;initial=[2.])
            attachment=joinpath(dir,"attachment.tsv")
            ROOTPlots.write_hist_draw_spec(attachment,h,spec,
                ROOTPlots.raw_hist_draw_limits(h,spec,GraphStyle(),"Hist"),
                "FUNC","Hist",GraphStyle();attached_fit=fit)
            write(attachment,replace(read(attachment,String),"LENTRY\tseries\t1"=>"LENTRY\tseries\t2"))
            failed=run(pipeline(ignorestatus(`$exe $attachment $(joinpath(dir,"bad-attachment"))`);
                stdout=devnull,stderr=devnull))
            @test failed.exitcode==2 # A controlled validation error, never SIGSEGV during TH1 cleanup.
            @test_throws ArgumentError render_root(h,fit,spec,joinpath(dir,"not-two-pad"))
        end
    end
end

@testset "predefined edge fits recover independent asymmetric fixtures" begin
    radial(x,p)=p[1]+p[2]/2*erfc((x-p[3])/(sqrt(2)*p[4]))
    double(x,p)=p[1]+p[2]/2*(erf((x-p[3])/(sqrt(2)*p[5]))-erf((x-p[4])/(sqrt(2)*p[6])))
    cases=((:erfc_edge,[1.7,24.,1.3,.18],collect(0.:.025:2.5),radial,
            [2.,22.,1.2,.23],[nothing,nothing,nothing,(.03,1.5)]),
           (:double_erf,[2.3,31.,-1.4,2.2,.21,.53],collect(-3.:.075:3.975),double,
            [2.,29.,-1.3,2.1,.3,.4],[nothing,nothing,(-2.,-1.),(1.5,2.8),(.03,1.5),(.03,1.5)]))
    if TEST_NATIVE_ROOT
        mktempdir() do dir
            for (model,truth,x,evaluate,initial,bounds) in cases
                y=evaluate.(x,Ref(truth)); errors=fill(.2,length(x))
                graph=GraphSeries(x,y;yerrors=errors)
                step=x[2]-x[1]
                edges=[x.-step/2; last(x)+step/2]
                for data in (graph,Hist1D(edges,y,errors.^2))
                    fitted=fit_root(data,model;initial,bounds)
                    @test fitted.status==0 && fitted.valid
                    @test fitted.parameters ≈ truth atol=.004
                    @test size(fitted.covariance)==(length(truth),length(truth))
                    curve=fit_curve(fitted;n=101)
                    @test curve.y ≈ evaluate.(curve.x,Ref(fitted.parameters)) rtol=1e-12
                end
                @test_throws ArgumentError fit_root(graph,model;initial)
                negative=copy(initial);negative[end]=0.
                @test_throws ArgumentError fit_root(graph,model;initial=negative,bounds)
                invalid_bounds=copy(bounds);invalid_bounds[end]=(0.,1.)
                @test_throws ArgumentError fit_root(graph,model;initial,bounds=invalid_bounds)
                width_indices=collect(ROOTPlots.edge_width_indices(model))
                fitted=fit_root(graph,model;initial=truth,fixed=width_indices)
                @test fitted.status==0 && fitted.parameters[width_indices]==truth[width_indices]
                fitted.parameters[end]=0.
                @test_throws ArgumentError fit_curve(fitted)
                @test_throws ArgumentError render_root(Hist1D(edges,y,errors.^2),GraphSpec(),joinpath(dir,"invalid");
                    draw="FUNC",attached_fit=fitted)
            end
            # A native caller must not bypass the Julia positive-width policy.
            x=collect(0.:.1:2.);g=GraphSeries(x,radial.(x,Ref(cases[1][2]));yerrors=fill(.2,length(x)))
            io=IOBuffer();println(io,"ROOTPlotsGraph\t1\nFIT\tgraph\terfc_edge\tchi2\teffective\t0\t2\t4")
            for (i,p) in enumerate(cases[1][2])
                println(io,"PARAM\t$i\t$p\t-Inf\tInf\t0")
            end
            ROOTPlots.write_graph_series(io,[g])
            input=joinpath(dir,"bad-width.tsv");write(input,String(take!(io)))
            exe=ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
            @test !success(pipeline(`$exe $input $(joinpath(dir,"bad-fit"))`;stdout=devnull,stderr=devnull))
        end
    end
end
