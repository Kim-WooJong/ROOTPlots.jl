abstract type SpecialGraph end

"""ROOT TGraphPolar. Angles and their errors are radians; radii must be nonnegative."""
struct PolarGraph <: SpecialGraph
    angle::Vector{Float64}
    radius::Vector{Float64}
    angle_errors::Vector{Float64}
    radius_errors::Vector{Float64}
    style::GraphStyle
end
PolarGraph(angle,radius;angle_errors=zeros(length(angle)),radius_errors=zeros(length(radius)),style=GraphStyle()) =
    validate_special(PolarGraph(Float64.(angle),Float64.(radius),Float64.(angle_errors),Float64.(radius_errors),style))

"""Native ROOT QQ plot: first sample on x, second on y (not paired observations)."""
struct QQGraph <: SpecialGraph
    x::Vector{Float64}
    y::Vector{Float64}
    style::GraphStyle
end
QQGraph(x,y;style=GraphStyle())=validate_special(QQGraph(Float64.(x),Float64.(y),style))

"""Native TSpline3/TSpline5 through ordered knots. No Julia interpolation approximation."""
struct SplineCurve <: SpecialGraph
    x::Vector{Float64}
    y::Vector{Float64}
    order::Int
    npx::Int
    style::GraphStyle
end
SplineCurve(x,y;order=3,npx=300,style=GraphStyle(color="#cc3311")) =
    validate_special(SplineCurve(Float64.(x),Float64.(y),Int(order),Int(npx),style))

"""Native ROOT TPie with explicit labels and colors."""
struct PieChart <: SpecialGraph
    values::Vector{Float64}
    labels::Vector{String}
    colors::Vector{String}
    style::GraphStyle
end
PieChart(values;labels=string.(1:length(values)),colors=fill("#4477aa",length(values)),style=GraphStyle()) =
    validate_special(PieChart(Float64.(values),String.(labels),String.(colors),style))

"""Native TEfficiency with integer counts and Clopper–Pearson confidence intervals."""
struct EfficiencyGraph <: SpecialGraph
    edges::Vector{Float64}
    passed::Vector{Float64}
    total::Vector{Float64}
    confidence::Float64
    style::GraphStyle
end
EfficiencyGraph(edges,passed,total;confidence=.682689492137,style=GraphStyle()) =
    validate_special(EfficiencyGraph(Float64.(edges),Float64.(passed),Float64.(total),Float64(confidence),style))

function validate_special(g::SpecialGraph)
    check_graph_style(g.style)
    if g isa PolarGraph
        n=length(g.angle)
        2<=n<=1_000_000 && all(length(a)==n for a in (g.radius,g.angle_errors,g.radius_errors)) &&
            all(isfinite,g.angle) && all(isfinite,g.radius) && all(isfinite,g.angle_errors) && all(isfinite,g.radius_errors) &&
            all(>=(0),g.radius) && all(>=(0),g.angle_errors) && all(>=(0),g.radius_errors) || throw(ArgumentError("invalid polar points/errors"))
    elseif g isa Union{QQGraph,SplineCurve}
        minimum_points=g isa SplineCurve ? (g.order==5 ? 5 : 3) : 2
        minimum_points<=length(g.x)<=1_000_000 && length(g.x)==length(g.y) && all(isfinite,g.x) && all(isfinite,g.y) || throw(ArgumentError("invalid graph samples"))
        if g isa SplineCurve
            g.order in (3,5) && all(diff(g.x).>0) && 2<=g.npx<=1_000_000 || throw(ArgumentError("spline needs ordered distinct x, order 3/5 and valid npx"))
        end
    elseif g isa PieChart
        1<=length(g.values)<=1000 && length(g.values)==length(g.labels)==length(g.colors) && all(isfinite,g.values) && all(>=(0),g.values) &&
            isfinite(sum(g.values)) && sum(g.values)>0 || throw(ArgumentError("pie needs matching nonnegative values with positive total"))
        foreach(scene_text,g.labels); foreach(scene_color,g.colors)
    else
        n=length(g.passed)
        1<=n<=1_000_000 && length(g.edges)==n+1 && length(g.total)==n && all(isfinite,g.edges) && all(diff(g.edges).>0) &&
            all(isfinite,g.passed) && all(isfinite,g.total) &&
            all(i->0<=g.passed[i]<=g.total[i]<=2.0^53 && isinteger(g.passed[i]) && isinteger(g.total[i]),1:n) &&
            isfinite(g.confidence) && 0<g.confidence<1 || throw(ArgumentError("efficiency needs ordered bins and integer 0<=passed<=total"))
    end
    return g
end

function write_special_plot(path,g::SpecialGraph,spec::GraphSpec;points=false)
    validate_special(g); validate_scene(spec;nseries=1)
    spec.legend in (:none,:topleft,:topright,:bottomleft,:bottomright) || throw(ArgumentError("invalid legend"))
    all(n->1<=n<=10000,spec.size) || throw(ArgumentError("invalid canvas size"))
    foreach(scene_text,(spec.title,spec.xlabel,spec.ylabel))
    polar=g isa PolarGraph; pie=g isa PieChart
    if polar || pie
        (spec.logx || spec.logy || spec.xlimits!==nothing || spec.ylimits!==nothing) &&
            throw(ArgumentError("polar/pie plots own non-Cartesian axes; Cartesian limits/logs are not applicable"))
    end
    kind=polar ? "polar" : pie ? "pie" : g isa QQGraph ? "qq" : g isa SplineCurve ? "spline" : "efficiency"
    n=polar ? length(g.angle) : pie ? length(g.values) : g isa EfficiencyGraph ? length(g.passed) : length(g.x)
    xrange=polar || pie ? (0.,1.) : g isa EfficiencyGraph ? (first(g.edges),last(g.edges)) : (minimum(g.x),maximum(g.x))
    yrange=polar || pie ? (0.,1.) : g isa EfficiencyGraph ? (0.,1.05) : (minimum(g.y),maximum(g.y))
    xrange[1]==xrange[2] && (xrange=(xrange[1]-.5,xrange[2]+.5))
    yrange[1]==yrange[2] && (yrange=(yrange[1]-.5,yrange[2]+.5))
    xr=spec.xlimits===nothing ? xrange : spec.xlimits
    yr=spec.ylimits===nothing ? (yrange[1]-.05*(yrange[2]-yrange[1]),yrange[2]+.05*(yrange[2]-yrange[1])) : spec.ylimits
    for (range,logscale) in ((xr,spec.logx),(yr,spec.logy))
        all(isfinite,range) && range[1]<range[2] && (!logscale || range[1]>0) || throw(ArgumentError("invalid special graph range"))
    end
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        println(io,join((points ? "SPECIALPOINTS" : "DRAWSPECIAL",spec.size...,kind,n,Int(spec.logx),Int(spec.logy),xr...,yr...,spec.legend,
                         g isa EfficiencyGraph ? g.confidence : .682689492137,g isa SplineCurve ? g.order : 3,g isa SplineCurve ? g.npx : 300),'\t'))
        println(io,join(("LABELS",scene_text(spec.title),scene_text(spec.xlabel),scene_text(spec.ylabel)),'\t'))
        write_scene(io,spec); panel_hist_style(io,g.style)
        for i in 1:n
            row=polar ? ("POINT",i,g.angle[i],g.radius[i],g.angle_errors[i],g.radius_errors[i]) :
                pie ? ("SLICE",i,g.values[i],scene_text(g.labels[i]),g.colors[i]) :
                g isa EfficiencyGraph ? ("BIN",i,g.edges[i],g.edges[i+1],g.passed[i],g.total[i]) : ("POINT",i,g.x[i],g.y[i])
            println(io,join(row,'\t'))
        end
        println(io,"END\t1")
    end
    return path
end

function render_root(g::SpecialGraph,spec::GraphSpec,stem::AbstractString;format=:png,root_config=nothing)
    return render_root(PlotPanel(g,spec),stem;format,root_config)
end

function special_points(g::Union{EfficiencyGraph,SplineCurve,QQGraph};root_config=nothing)
    validate_special(g)
    config=find_root_config(root_config); config===nothing && throw(ArgumentError("local native ROOT not found"))
    mktempdir() do dir
        input=write_special_plot(joinpath(dir,"points.tsv"),g,GraphSpec();points=true)
        exe=graphfit_executable(config); output=joinpath(dir,"result.tsv")
        run(Cmd([exe,input,output]))
        rows=readlines(output)
        first(rows)=="ROOTPlotsSpecialPoints\t1" && last(rows)=="END" || error("invalid special point result")
        samples=[parse.(Float64,split(row,'\t')[2:end]) for row in rows[2:end-1]]
        xx=[r[1] for r in samples]; yy=[r[2] for r in samples]
        g isa Union{SplineCurve,QQGraph} && return GraphSeries(xx,yy;draw=g isa SplineCurve ? :line : :points,errors=:none,style=g.style)
        return GraphSeries(xx,yy;xerrors=(low=[r[3] for r in samples],high=[r[4] for r in samples]),
            yerrors=(low=[r[5] for r in samples],high=[r[6] for r in samples]),style=g.style)
    end
end
spline_points(g::SplineCurve;kwargs...)=special_points(g;kwargs...)
efficiency_points(g::EfficiencyGraph;kwargs...)=special_points(g;kwargs...)
qq_points(g::QQGraph;kwargs...)=special_points(g;kwargs...)
