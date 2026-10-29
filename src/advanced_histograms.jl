"""Owned THStack inputs. True stacks require identical bins and nonnegative contents.

Use `draw="NOSTACK HIST"` for ROOT's unstacked THStack view. Fill attributes
are separate from GraphStyle so existing graph styles and protocols stay unchanged.
"""
struct HistStack
    hists::Vector{Hist1D}
    labels::Vector{String}
    styles::Vector{GraphStyle}
    fill_colors::Vector{String}
    fill_styles::Vector{Int}
    draw::String
    function HistStack(hists::Vector{Hist1D},labels::Vector{String},styles::Vector{GraphStyle},
                       fill_colors::Vector{String},fill_styles::Vector{Int},draw::String)
        s=new([Hist1D(h.edges,h.sumw,h.sumw2) for h in hists],copy(labels),copy(styles),
              copy(fill_colors),copy(fill_styles),draw)
        validate_stack(s)
        return s
    end
end

function HistStack(hists::AbstractVector{Hist1D}; labels=fill("",length(hists)),
                   styles=fill(GraphStyle(),length(hists)),
                   fill_colors=[s.color for s in styles],
                   fill_styles=fill(1001,length(hists)), draw="HIST")
    stack=HistStack(collect(hists),String.(labels),
                    collect(styles),String.(fill_colors),Int.(fill_styles),String(draw))
    validate_stack(stack)
    return stack
end

function validate_stack(s::HistStack;draw=s.draw)
    n=length(s.hists)
    1 <= n <= 64 && length(s.labels)==length(s.styles)==length(s.fill_colors)==length(s.fill_styles)==n ||
        throw(ArgumentError("stack needs 1:64 matching histograms, labels, and styles"))
    root_draw_option(draw)
    unstacked=occursin("NOSTACK",uppercase(draw))
    for i in 1:n
        validate_stats_hist(s.hists[i]); check_graph_style(s.styles[i])
        length(s.hists[i].sumw)<=1_000_000 || throw(ArgumentError("too many stack bins"))
        scene_text(s.labels[i]); scene_color(s.fill_colors[i])
        s.fill_styles[i] in (0,1001) || 3000 <= s.fill_styles[i] <= 3999 ||
            throw(ArgumentError("fill style must be 0, 1001, or 3000:3999"))
        if !unstacked
            s.hists[i].edges==s.hists[1].edges || throw(ArgumentError("stacked histograms need identical edges"))
            all(>=(0),s.hists[i].sumw) || throw(ArgumentError("negative bins need NOSTACK"))
        end
    end
    return s
end

abstract type AbstractRootProfile end

"""Owned samples for a genuine ROOT TProfile. Nonnegative weights; visible bins are half-open.

`error` is :mean (standard error), :spread, :integer, or :gaussian (1/√Σw).
Under/overflow samples are retained for ROOT filling but omitted by profile_stats.
"""
struct Profile1D <: AbstractRootProfile
    xedges::Vector{Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    weights::Vector{Float64}
    error::Symbol
    function Profile1D(edges::AbstractVector,x::AbstractVector,y::AbstractVector;
                       weights=nothing,error=:mean)
        p=new(Float64.(edges),Float64.(x),Float64.(y),
              weights===nothing ? ones(length(x)) : Float64.(weights),error)
        validate_profile(p)
        return p
    end
end

"""Owned samples for native TProfile2D: average z at each (x,y) bin."""
struct Profile2D <: AbstractRootProfile
    xedges::Vector{Float64}
    yedges::Vector{Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    weights::Vector{Float64}
    error::Symbol
    function Profile2D(xedges::AbstractVector,yedges::AbstractVector,
                       x::AbstractVector,y::AbstractVector,z::AbstractVector;
                       weights=nothing,error=:mean)
        p=new(Float64.(xedges),Float64.(yedges),Float64.(x),Float64.(y),Float64.(z),
              weights===nothing ? ones(length(x)) : Float64.(weights),error)
        validate_profile(p)
        return p
    end
end

"""Owned samples for native TProfile3D: average t at each (x,y,z) bin."""
struct Profile3D <: AbstractRootProfile
    xedges::Vector{Float64}
    yedges::Vector{Float64}
    zedges::Vector{Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    t::Vector{Float64}
    weights::Vector{Float64}
    error::Symbol
    function Profile3D(xedges::AbstractVector,yedges::AbstractVector,zedges::AbstractVector,
                       x::AbstractVector,y::AbstractVector,z::AbstractVector,t::AbstractVector;
                       weights=nothing,error=:mean)
        p=new(Float64.(xedges),Float64.(yedges),Float64.(zedges),Float64.(x),Float64.(y),
              Float64.(z),Float64.(t),weights===nothing ? ones(length(x)) : Float64.(weights),error)
        validate_profile(p)
        return p
    end
end

profile_edges(p::Profile1D)=(p.xedges,)
profile_edges(p::Profile2D)=(p.xedges,p.yedges)
profile_edges(p::Profile3D)=(p.xedges,p.yedges,p.zedges)
profile_columns(p::Profile1D)=(p.x,p.y)
profile_columns(p::Profile2D)=(p.x,p.y,p.z)
profile_columns(p::Profile3D)=(p.x,p.y,p.z,p.t)

function validate_profile(p::AbstractRootProfile)
    p.error in (:mean,:spread,:integer,:gaussian) || throw(ArgumentError("unknown profile error mode"))
    edges=profile_edges(p); columns=profile_columns(p)
    n=length(p.weights)
    0 <= n <= 1_000_000 && all(length(c)==n for c in columns) ||
        throw(ArgumentError("profile needs matching arrays with at most 1,000,000 samples"))
    cells=1
    for e in edges
        cap=length(edges)==3 ? 256 : 2048
        2 <= length(e) <= cap+1 && all(isfinite,e) && all(diff(e).>0) ||
            throw(ArgumentError("profile edges must be finite, increasing, and within axis limits"))
        cells*=length(e)-1
    end
    cells<=1_000_000 || throw(ArgumentError("profile has too many cells"))
    all(isfinite,p.weights) && all(>=(0),p.weights) && all(c->all(isfinite,c),columns) ||
        throw(ArgumentError("profile samples must be finite with nonnegative weights"))
    p.error==:gaussian && !all(>(0),p.weights) && throw(ArgumentError("gaussian profile errors need positive weights"))
    for i in 1:n
        w=p.weights[i]
        isfinite(w*w) && all(c->isfinite(w*c[i]) && isfinite(w*c[i]*c[i]),columns) ||
            throw(ArgumentError("profile weighted moments overflow Float64"))
    end
    return p
end

"""Owned pre-binned TH2Poly cells. Each polygon is an N×2 matrix or (x,y) vectors.

Any simple polygon (convex or non-convex, e.g. L-shaped or star-shaped) without
holes is accepted; see `simple_polygon` for the vertex checks (3:256 finite
vertices, non-zero area, no self-intersection; a repeated closing vertex is
dropped). Bins may overlap; ROOT paints them in input order. No
event-filling/overlap assignment is performed.
"""
struct PolygonHist2D
    polygons::Vector{Matrix{Float64}}
    values::Vector{Float64}
    function PolygonHist2D(polygons::AbstractVector,values::AbstractVector)
        owned=Matrix{Float64}[]
        for polygon in polygons
            if polygon isa AbstractMatrix
                size(polygon,2)==2 || throw(ArgumentError("polygon matrix must be N×2"))
                xx,yy=simple_polygon(polygon[:,1],polygon[:,2])
                push!(owned,hcat(xx,yy))
            elseif polygon isa Tuple && length(polygon)==2
                xx,yy=simple_polygon(polygon[1],polygon[2])
                push!(owned,hcat(xx,yy))
            else
                throw(ArgumentError("polygon must be an N×2 matrix or (x,y) vectors"))
            end
        end
        p=new(owned,Float64.(values))
        validate_poly(p)
        return p
    end
end

function validate_poly(p::PolygonHist2D)
    1 <= length(p.polygons) <= 10_000 && length(p.polygons)==length(p.values) && all(isfinite,p.values) ||
        throw(ArgumentError("polygon histogram needs 1:10000 finite cell values"))
    sum(size(v,1) for v in p.polygons)<=100_000 || throw(ArgumentError("too many polygon vertices"))
    for v in p.polygons
        size(v,2)==2 || throw(ArgumentError("polygon needs N×2 xy vertices"))
        n=size(v,1)
        # Stored polygons already had any closing duplicate removed; a changed count means mutation.
        length(first(simple_polygon(v[:,1],v[:,2])))==n || throw(ArgumentError("polygon vertices changed"))
    end
    return p
end

function advanced_header(io,command,data,spec;draw,stats=false)
    dim=data isa AbstractRootProfile ? length(profile_edges(data)) : data isa HistStack ? 1 : 2
    logz=spec isa GraphSpec ? false : spec.logz
    palette=spec isa GraphSpec ? :viridis : spec.palette
    contours=spec isa Graph2DSpec ? spec.contours : 10
    theta=spec isa Graph3DSpec ? spec.theta : 30.
    phi=spec isa Graph3DSpec ? spec.phi : 30.
    xl=spec isa GraphSpec && spec.xlimits!==nothing ? spec.xlimits : (0.,1.)
    yl=spec isa GraphSpec && spec.ylimits!==nothing ? spec.ylimits : (0.,1.)
    count=data isa HistStack ? length(data.hists) : data isa PolygonHist2D ? length(data.polygons) : length(data.weights)
    error=data isa AbstractRootProfile ? String(data.error) : "-"
    println(io,join((command,spec.size...,dim,root_draw_option(draw),Int(spec.logx),Int(spec.logy),Int(logz),
        palette,contours,theta,phi,Int(spec isa GraphSpec && spec.xlimits!==nothing),xl...,
        Int(spec isa GraphSpec && spec.ylimits!==nothing),yl...,Int(stats),count,error,
        spec isa GraphSpec ? String(spec.legend) : "none"),'\t'))
    zlabel=spec isa GraphSpec ? "" : spec.zlabel
    println(io,join(("LABELS",(scene_text(s) for s in (spec.title,spec.xlabel,spec.ylabel,zlabel))...),'\t'))
    spec isa GraphSpec && write_scene(io,spec)
end

function write_advanced_spec(path,data,spec;draw,stats=false,command=nothing)
    command===nothing && (command=data isa HistStack ? "DRAWSTACK" : data isa PolygonHist2D ? "DRAWPOLY" : "DRAWPROFILE")
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        advanced_header(io,command,data,spec;draw,stats)
        if data isa HistStack
            for (k,h) in enumerate(data.hists)
                s=data.styles[k]
                println(io,join(("HLAYER",k,scene_text(data.labels[k]),s.color,s.marker_style,s.marker_size,
                    s.line_width,s.line_style,data.fill_colors[k],data.fill_styles[k]),'\t'))
                println(io,"HIST\t",length(h.sumw))
                for (i,e) in enumerate(h.edges) println(io,"EDGE\t",i,'\t',e) end
                for i in eachindex(h.sumw) println(io,"BIN\t",i,'\t',h.sumw[i],'\t',h.sumw2[i]) end
            end
        elseif data isa AbstractRootProfile
            for (axis,edges) in enumerate(profile_edges(data))
                println(io,"PEDGES\t",axis,'\t',length(edges))
                for (i,e) in enumerate(edges) println(io,"PEDGE\t",axis,'\t',i,'\t',e) end
            end
            for i in eachindex(data.weights)
                println(io,join(("SAMPLE",i,(c[i] for c in profile_columns(data))...,data.weights[i]),'\t'))
            end
        else
            for (i,v) in enumerate(data.polygons)
                println(io,"POLYGON\t",i,'\t',size(v,1),'\t',data.values[i])
                for j in axes(v,1) println(io,"VERTEX\t",j,'\t',v[j,1],'\t',v[j,2]) end
            end
        end
        count=data isa HistStack ? length(data.hists) : data isa PolygonHist2D ? length(data.polygons) : length(data.weights)
        println(io,"END\t",count)
    end
    return path
end

function validate_advanced(data,spec,draw)
    root_draw_option(draw)
    compatible=(data isa Union{HistStack,Profile1D} && spec isa GraphSpec) ||
        (data isa Profile2D && spec isa Union{Graph2DSpec,Graph3DSpec}) ||
        (data isa Profile3D && spec isa Graph3DSpec) ||
        (data isa PolygonHist2D && spec isa Graph2DSpec)
    compatible || throw(ArgumentError("advanced histogram and plot specification dimensions do not match"))
    if spec isa GraphSpec
        validate_scene(spec;nseries=data isa HistStack ? length(data.hists) : 1)
        spec.legend in (:topright,:topleft,:bottomright,:bottomleft,:none) || throw(ArgumentError("invalid legend"))
        all(n->1<=n<=10000,spec.size) || throw(ArgumentError("invalid canvas size"))
        for (limits,logscale) in ((spec.xlimits,spec.logx),(spec.ylimits,spec.logy))
            limits===nothing || (all(isfinite,limits) && limits[1]<limits[2] && (!logscale || limits[1]>0)) ||
                throw(ArgumentError("invalid axis limits"))
        end
    else
        validate_nd_spec(spec)
    end
    for s in (spec.title,spec.xlabel,spec.ylabel) scene_text(s) end
    data isa HistStack ? validate_stack(data;draw) : data isa PolygonHist2D ? validate_poly(data) : validate_profile(data)
    edges=data isa HistStack ? (minimum(first(h.edges) for h in data.hists),) :
        data isa PolygonHist2D ? (minimum(minimum(v[:,1]) for v in data.polygons),minimum(minimum(v[:,2]) for v in data.polygons)) :
        first.(profile_edges(data))
    logs=spec isa GraphSpec ? (spec.logx,) : (spec.logx,spec.logy,spec.logz)
    for i in eachindex(edges)
        logs[i] && edges[i]<=0 && throw(ArgumentError("log coordinate axis needs positive bin edges"))
    end
    return data
end

"""Validate and serialize a standalone advanced plot without executing ROOT (layout helper)."""
function write_advanced_plot(path,data,spec;draw=nothing,stats=false)
    option=if draw!==nothing
        String(draw)
    elseif data isa HistStack
        data.draw
    elseif data isa Profile1D
        "E1"
    elseif spec.draw isa AbstractString
        spec.draw
    elseif spec isa Graph2DSpec
        spec.draw==:heatmap ? "COLZ" : "CONT3"
    elseif data isa Profile2D && spec.draw==:surface
        "SURF2Z"
    else
        throw(ArgumentError("this advanced plot needs a ROOT draw option string"))
    end
    validate_advanced(data,spec,option)
    return write_advanced_spec(path,data,spec;draw=option,stats)
end

function render_advanced(data,spec,outstem;draw,format=:png,stats=false,root_config=nothing)
    paths=plot_output_paths(outstem,format)
    validate_advanced(data,spec,draw)
    config=find_root_config(root_config)
    config===nothing && throw(ArgumentError("local native ROOT not found"))
    exe=graphfit_executable(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        input=write_advanced_spec(joinpath(dir,"advanced.tsv"),data,spec;draw,stats)
        stem=joinpath(dir,"rendered")
        run(Cmd([exe,input,stem]))
        move_plot_outputs(stem,paths)
    end
    return paths
end

render_root(data::HistStack,spec::GraphSpec,outstem::AbstractString;
            draw=data.draw,format=:png,stats=false,root_config=nothing)=
    render_advanced(data,spec,outstem;draw=String(draw),format,stats,root_config)
render_root(data::Profile1D,spec::GraphSpec,outstem::AbstractString;
            draw="E1",format=:png,stats=false,root_config=nothing)=
    render_advanced(data,spec,outstem;draw=String(draw),format,stats,root_config)
function render_root(data::Union{Profile2D,PolygonHist2D},spec::Graph2DSpec,outstem::AbstractString;
                     format=:png,stats=false,root_config=nothing)
    draw=spec.draw isa AbstractString ? spec.draw : spec.draw==:heatmap ? "COLZ" : "CONT3"
    return render_advanced(data,spec,outstem;draw,format,stats,root_config)
end
function render_root(data::Profile2D,spec::Graph3DSpec,outstem::AbstractString;
                     format=:png,stats=false,root_config=nothing)
    spec.draw isa AbstractString || spec.draw==:surface || throw(ArgumentError("Profile2D needs :surface or ROOT draw string"))
    draw=spec.draw isa AbstractString ? spec.draw : "SURF2Z"
    return render_advanced(data,spec,outstem;draw,format,stats,root_config)
end
function render_root(data::Profile3D,spec::Graph3DSpec,outstem::AbstractString;
                     format=:png,stats=false,root_config=nothing)
    spec.draw isa AbstractString || throw(ArgumentError("Profile3D needs a ROOT draw string, e.g. BOX2Z"))
    return render_advanced(data,spec,outstem;draw=spec.draw,format,stats,root_config)
end

"""Read native visible-bin profile means/errors without a ROOT file.

Arrays use Julia axis order. `entries` is Σw, not the number of samples;
`effective_entries` is (Σw)²/Σw². Empty cells return zero. Native ROOT computes
all quantities, including the selected error mode.
"""
function profile_stats(p::AbstractRootProfile;root_config=nothing)
    validate_profile(p)
    config=find_root_config(root_config)
    config===nothing && throw(ArgumentError("local native ROOT not found"))
    exe=graphfit_executable(config)
    spec=p isa Profile1D ? GraphSpec(legend=:none) : p isa Profile2D ? Graph2DSpec() : Graph3DSpec(draw="BOX2Z")
    shape=Tuple(length(e)-1 for e in profile_edges(p))
    arrays=ntuple(_->zeros(shape),4)
    mktempdir() do dir
        input=write_advanced_spec(joinpath(dir,"profile.tsv"),p,spec;draw="",command="PROFILESTATS")
        output=joinpath(dir,"summary.tsv")
        run(Cmd([exe,input,output]))
        lines=readlines(output)
        length(lines)==prod(shape)+1 && first(lines)=="ROOTPlotsProfile\t1" || error("invalid native profile summary")
        for (i,line) in enumerate(Iterators.drop(lines,1))
            fields=split(line,'\t')
            length(fields)==6 && fields[1]=="BIN" && parse(Int,fields[2])==i || error("invalid profile summary row")
            for k in 1:4 arrays[k][i]=parse(Float64,fields[k+2]) end
        end
    end
    return (mean=arrays[1],error=arrays[2],entries=arrays[3],effective_entries=arrays[4])
end
