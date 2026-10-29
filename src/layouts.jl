"""
    PlotPanel(data, spec; primitives=PlotPrimitive[], kwargs...)

One native ROOT pad. Data/spec are read-shared; do not mutate them while rendering.

Axis/canvas keywords for histogram, grid and graph-layer panels: `xaxis`, `yaxis`,
`zaxis` ([`AxisStyle`](@ref)) and `canvas` ([`CanvasStyle`](@ref)) replace the spec's
fields of the same name; `xlimits`, `ylimits`, `zlimits` set the drawn ranges.

Statistics/fit-box keywords for histogram and graph panels: `stats` (box of the
first histogram layer; Hist1D, Hist2D, Hist3D, profiles, TH2Poly; see
[`stats_option`](@ref)), `fitstats` and `attached_fit` (a `FitResult` attached
as a native TF1 to the first `Hist1D` layer, or else the first `GraphSeries`
layer; see [`fit_option`](@ref)), `fit_label`, `stats_position`,
`stats_name` (default: the layer label), and `stats_headroom` (1-D panels:
raise the automatic y maximum so the data stay below the box). Each panel keeps
its own box settings.
"""
struct PlotPanel{D,S,K}
    data::D
    spec::S
    options::K
    primitives::Vector{PlotPrimitive}
end
PlotPanel(data, spec; primitives=PlotPrimitive[], kwargs...) =
    PlotPanel(data, spec, (; kwargs...), PlotPrimitive[primitives...])

"""
    PlotFigure(cells; layout=(1, length(cells)), size=(1000, 700), widths=nothing, heights=nothing,
               style=nothing)
    PlotFigure(cells::AbstractMatrix; size=(1000, 700), widths=nothing, heights=nothing, style=nothing)

Grid of native ROOT pads on one canvas, in row-major order (no raster
composition). Each cell is a [`PlotPanel`](@ref) or another `PlotFigure`; a
nested figure is drawn as a sub-grid inside its cell (ROOT `TPad::Divide` of that
cell's pad), so sub-grids can be nested up to 4 levels deep. The `size` of a
nested figure is ignored; only the outermost figure sets the canvas size.

`widths` (one positive weight per column) and `heights` (one per row) make
columns/rows proportionally wider/taller; omitted weights give ROOT's regular
`Divide` grid. A figure holds at most 64 panels in total, with 1:16 rows/columns
per grid.

`style` selects one ROOT `TStyle` for the whole canvas (`gROOT->SetStyle` is
process-global): it is given to every panel whose spec has no `style`, so with
`style=:CMS` each such panel also gets the CMS label (`cms_label=false` in a
panel spec switches it off). Panel specs that name a style must all agree.
"""
struct PlotFigure
    panels::Vector{Union{PlotPanel,PlotFigure}}
    layout::Tuple{Int,Int}
    size::Tuple{Int,Int}
    widths::Union{Nothing,Vector{Float64}}
    heights::Union{Nothing,Vector{Float64}}
    style::Union{Nothing,Symbol}
    function PlotFigure(panels; layout=(1,length(panels)), size=(1000,700), widths=nothing, heights=nothing,
                        style=nothing)
        all(p -> p isa Union{PlotPanel,PlotFigure}, panels) ||
            throw(ArgumentError("figure cells must be PlotPanel or PlotFigure values"))
        p = Union{PlotPanel,PlotFigure}[panels...]
        rows, cols = layout
        rows isa Integer && cols isa Integer && 1 <= rows <= 16 && 1 <= cols <= 16 &&
            1 <= length(p) <= 64 && length(p) <= rows*cols <= 64 ||
            throw(ArgumentError("layout needs 1:64 pads in a 1:16 row/column grid"))
        length(size) == 2 && all(v -> v isa Integer && 100 <= v <= 10000, size) ||
            throw(ArgumentError("figure size must be two integers in 100:10000"))
        weights(w, n, name) = w === nothing ? nothing :
            (length(w) == n && all(v -> v isa Real && isfinite(v) && v > 0, w) && maximum(w)/minimum(w) <= 100 ?
             Float64.(collect(w)) : throw(ArgumentError("$name needs $n finite positive weights (max/min ≤ 100)")))
        check_style_name(style)
        figure = new(p, (Int(rows),Int(cols)), (Int(size[1]),Int(size[2])),
                     weights(widths, Int(cols), "widths"), weights(heights, Int(rows), "heights"), style)
        figure_depth(figure) <= 4 || throw(ArgumentError("figures can be nested at most 4 levels deep"))
        length(figure_leaves(figure)) <= 64 || throw(ArgumentError("a figure holds at most 64 panels in total"))
        return figure
    end
end
PlotFigure(panels::AbstractMatrix; size=(1000,700), widths=nothing, heights=nothing, style=nothing) =
    PlotFigure(vec(permutedims(panels)); layout=Base.size(panels), size, widths, heights, style)

"""
Give the figure's `style` (inherited by nested figures) to every panel spec without
one, then check that all panels name the same style (gROOT->SetStyle is global).
"""
function styled_figure(f::PlotFigure, inherited=nothing)
    style = f.style === nothing ? inherited : f.style
    cells = map(f.panels) do cell
        cell isa PlotFigure && return styled_figure(cell, style)
        (style === nothing || !hasproperty(cell.spec, :style) || cell.spec.style !== nothing) && return cell
        PlotPanel(cell.data, with_spec_fields(cell.spec; style), cell.options, cell.primitives)
    end
    return PlotFigure(cells; layout=f.layout, size=f.size, widths=f.widths, heights=f.heights, style=f.style)
end
figure_styles(f::PlotFigure) = vcat(f.style === nothing ? Symbol[] : [f.style],
    (c isa PlotFigure ? figure_styles(c) :
     hasproperty(c.spec, :style) && c.spec.style !== nothing ? [c.spec.style] : Symbol[] for c in f.panels)...)
function figure_style(f::PlotFigure)
    names = unique(lowercase.(String.(figure_styles(f))))
    length(names) <= 1 ||
        throw(ArgumentError("panels request different ROOT styles ($(join(names, ", "))); gROOT->SetStyle is process-global"))
    return f
end

figure_depth(f::PlotFigure) = 1 + maximum((c isa PlotFigure ? figure_depth(c) : 0) for c in f.panels)
figure_leaves(f::PlotFigure) = reduce(vcat, (c isa PlotFigure ? figure_leaves(c) : PlotPanel[c] for c in f.panels);
                                      init=PlotPanel[])

function panel_scene(panel, nseries)
    spec, options = panel.spec, panel.options
    common = (:xaxis,:yaxis,:canvas,:lines,:texts,:legend_spec,:legend,:xlimits,:ylimits)
    if spec isa GraphSpec
        scene = with_graph_spec(spec; (key=>options[key] for key in common if haskey(options,key))...)
    else
        defaults = (xlabel=spec.xlabel,ylabel=spec.ylabel,title=spec.title,
                    logx=spec.logx,logy=spec.logy,size=spec.size,legend=:none,
                    xaxis=spec.xaxis,yaxis=spec.yaxis,
                    style=spec.style,nd_label_fields(spec)...,
                    # The native painter picks the right margin (0.18, or 0.2 for a titled palette,
                    # widened until the palette labels, exponent and z title fit);
                    # explicit Graph2DSpec/Graph3DSpec canvas fields override the layout defaults.
                    # A selected ROOT style keeps its own left/bottom margins.
                    canvas=spec.style === nothing ?
                           merge_canvas(CanvasStyle(left_margin=0.14,bottom_margin=0.13),spec.canvas) : spec.canvas)
        overrides = (; (key=>options[key] for key in common if haskey(options,key))...)
        scene = GraphSpec(; merge(defaults, overrides)...)
    end
    validate_scene(scene; nseries)
    return scene
end

panel_layers(data::Union{Hist1D,GraphSeries,BentErrorSeries,MultiErrorSeries,Hist2D,Hist3D,Grid2D,Grid3D}) = [data]
function panel_layers(data::AbstractVector)
    all(x->x isa Union{Hist1D,GraphSeries,ModelCurve,BentErrorSeries,MultiErrorSeries},data) ||
        throw(ArgumentError("XY overlay layers must be Hist1D, GraphSeries, BentErrorSeries, MultiErrorSeries or ModelCurve"))
    return collect(data)
end
function panel_layers(data::Tuple{Hist1D,<:AbstractVector{GraphSeries}})
    return [data[1]; data[2]]
end

function panel_limits(layers, scene)
    xlo, xhi, ylo, yhi = Inf, -Inf, Inf, -Inf
    positive = Float64[]
    for layer in layers
        if layer isa Hist1D
            validate_stats_hist(layer)
            xlo=min(xlo,first(layer.edges)); xhi=max(xhi,last(layer.edges))
            low=layer.sumw.-sqrt.(layer.sumw2); high=layer.sumw.+sqrt.(layer.sumw2)
            ylo=min(ylo,minimum(low),0.); yhi=max(yhi,maximum(high),0.)
            append!(positive,filter(>(0),vcat(layer.sumw,low,high)))
        else
            range_scene=with_graph_spec(scene;legend_spec=nothing,lines=PlotLine[],texts=PlotText[])
            limits=layer isa ModelCurve ? model_curve_limits(layer,range_scene) :
                   layer isa ErrorGraph ? error_graph_limits(layer,range_scene) : graph_limits([layer],range_scene)
            xlo=min(xlo,limits[1][1]); xhi=max(xhi,limits[1][2])
            ylo=min(ylo,limits[2][1]); yhi=max(yhi,limits[2][2])
            layer isa Union{GraphSeries,ErrorGraph} && append!(positive,filter(>(0),layer.y))
            layer isa ModelCurve && append!(positive,filter(>(0),[limits[2]...]))
        end
    end
    xrange=scene.xlimits === nothing ? (xlo,xhi) : scene.xlimits
    yrange = if scene.ylimits !== nothing
        scene.ylimits
    elseif scene.logy
        isempty(positive) && throw(ArgumentError("log y needs positive visible values"))
        (minimum(positive)/2,max(yhi,maximum(positive))*1.15)
    else
        (ylo,yhi+max(1.,yhi-ylo)*0.12)
    end
    for (limits,logscale) in ((xrange,scene.logx),(yrange,scene.logy))
        all(isfinite,limits) && limits[1]<limits[2] && (!logscale || limits[1]>0) ||
            throw(ArgumentError("invalid panel limits"))
    end
    return (xrange,yrange)
end

function panel_hist_style(io, style)
    check_graph_style(style)
    println(io,join(("STYLE",style.color,style.marker_style,style.marker_size,
                    style.line_width,style.line_style,
                    something(style.fill_color,style.color),style.fill_style,style.fill_alpha),'\t'))
end

function write_basic_panel(path,panel)
    layers=panel_layers(panel.data)
    1 <= length(layers) <= 64 || throw(ArgumentError("panel needs 1:64 layers"))
    nd=any(x -> x isa Union{Hist2D,Hist3D,Grid2D,Grid3D},layers)
    nd && length(layers)!=1 && throw(ArgumentError("ND layers cannot be mixed in this layout API"))
    spec,opts=panel.spec,panel.options
    allowed=(:draw,:style,:styles,:label,:labels,:xaxis,:yaxis,:zaxis,:canvas,:lines,:texts,
             :legend_spec,:legend,:xlimits,:ylimits,:zlimits,:contour_levels,
             :stats,:fitstats,:stats_position,:stats_name,:stats_headroom,:attached_fit,:fit_label)
    unknown=filter(k -> !(k in allowed),collect(keys(opts)))
    isempty(unknown) || throw(ArgumentError("unknown panel keyword: $(join(unknown, ", "))"))
    scene=panel_scene(panel,length(layers))
    zaxis=get(opts,:zaxis,hasproperty(spec,:zaxis) ? spec.zaxis : AxisStyle())
    zaxis isa AxisStyle || throw(ArgumentError("zaxis must be an AxisStyle"))
    validate_axis(zaxis; name=:z)
    levels=Float64.(get(opts,:contour_levels,Float64[]))
    length(levels)<=255 && all(isfinite,levels) && all(diff(levels).>0) ||
        throw(ArgumentError("contour levels must be finite and strictly increasing (at most 255)"))
    palette_name=nd ? spec.palette : :viridis
    theta=spec isa Graph3DSpec ? spec.theta : 30.
    phi=spec isa Graph3DSpec ? spec.phi : 30.
    logz=nd && spec.logz
    counts=nd ? spec.contours : 10
    nd && validate_nd_spec(spec)
    dimensions=nd ? (first(layers) isa Union{Hist3D,Grid3D} ? 3 : 2) : 1
    ranges = if nd
        h=first(layers)
        h isa Union{Hist2D,Hist3D} ? validate_stats_hist_nd(h) : validate_grid(h)
        xr=scene.xlimits === nothing ? (first(h.xedges),last(h.xedges)) : scene.xlimits
        yr=scene.ylimits === nothing ? (first(h.yedges),last(h.yedges)) : scene.ylimits
        zr=get(opts,:zlimits,nothing)
        zdefault=dimensions==3 ? (first(h.zedges),last(h.zedges)) : (0.,1.)
        (xr,yr,zr === nothing ? zdefault : Tuple(Float64.(zr)))
    else
        xy=panel_limits(layers,scene)
        (xy[1],xy[2],(0.,1.))
    end
    # Bin labels need the histogram bins of the drawn frame: the first Hist1D layer of a
    # 1-D panel (its frame takes that binning) or the 2-D/3-D histogram itself.
    let h=first(layers)
        bins = if nd
            (length(h.xedges)-1, length(h.yedges)-1, dimensions==3 ? length(h.zedges)-1 : nothing)
        else
            k=findfirst(x -> x isa Hist1D, layers)
            (k === nothing ? nothing : length(layers[k].sumw), nothing, nothing)
        end
        what = nd ? "this panel" : "a panel without a Hist1D layer"
        check_axis_bins(scene.xaxis,bins[1],:x,what); check_axis_bins(scene.yaxis,bins[2],:y,what)
        check_axis_bins(zaxis,bins[3],:z,what)
    end
    for (axis, limits, logscale) in ((:x,ranges[1],scene.logx),(:y,ranges[2],scene.logy),(:z,ranges[3],logz))
        # For 2D histograms z is the color value, not a spatial bin axis.
        (axis==:z && dimensions==2 && !haskey(opts,:zlimits)) && continue
        length(limits)==2 && all(isfinite,limits) && limits[1]<limits[2] && (!logscale || limits[1]>0) ||
            throw(ArgumentError("invalid $axis panel limits"))
    end
    labels=get(opts,:labels,[get(opts,:label,h isa Union{GraphSeries,ModelCurve,ErrorGraph} ? h.label : "") for h in layers])
    styles=get(opts,:styles,[get(opts,:style,h isa Union{GraphSeries,ModelCurve,ErrorGraph} ? h.style : GraphStyle()) for h in layers])
    length(labels)==length(styles)==length(layers) || throw(ArgumentError("layer labels/styles must match"))
    # Statistics box of the first histogram layer; fit attached to the first Hist1D, else GraphSeries.
    first_hist=findfirst(h -> h isa Union{Hist1D,Hist2D,Hist3D},layers)
    attached_fit=get(opts,:attached_fit,nothing)
    fit_target=nothing
    if attached_fit !== nothing
        validate_attached_fit(attached_fit)
        fit_target=something(findfirst(h -> h isa Hist1D,layers),findfirst(h -> h isa GraphSeries,layers),Some(nothing))
        fit_target === nothing && throw(ArgumentError("attached_fit needs a Hist1D or GraphSeries layer"))
    end
    owner=first_hist !== nothing ? first_hist : fit_target
    request=stats_request(; stats=get(opts,:stats,false),fitstats=get(opts,:fitstats,false),
                          stats_position=get(opts,:stats_position,:auto),fit=attached_fit,
                          stats_headroom=get(opts,:stats_headroom,false),ylimits=scene.ylimits,
                          one_dimensional=!nd,
                          allow_stats=first_hist !== nothing,what="panels without a histogram layer",
                          name=stats_header(get(opts,:stats_name,nothing),owner === nothing ? "" : String(labels[owner])))
    fit_label=String(get(opts,:fit_label,""))
    clean_field(fit_label)
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        println(io,join(("PLOT",spec.size...,dimensions,Int(scene.logx),Int(scene.logy),Int(logz),
                         theta,phi,palette_name,counts,length(layers),scene.legend),'\t'))
        write_stats_row(io,request)
        if attached_fit !== nothing
            println(io,join(("FITTARGET",fit_target,scene_text(fit_label)),'\t'))
            write_attached_fit(io,attached_fit)
        end
        println(io,join(("LABELS",(scene_text(s) for s in (scene.title,scene.xlabel,scene.ylabel,nd ? spec.zlabel : ""))...),'\t'))
        println(io,join(("LIMITS",1,ranges[1]...,1,ranges[2]...,Int(haskey(opts,:zlimits)),ranges[3]...),'\t'))
        write_scene(io,scene)
        write_axis(io,"ZSTYLE","z",zaxis)
        println(io,join(("LEVELS",length(levels),levels...),'\t'))
        for (i,h) in enumerate(layers)
            label=scene_text(String(labels[i]))
            if h isa ModelCurve
                validate_model_curve(h); check_graph_style(styles[i])
                println(io,"LAYER\tFUNCTION\t",label)
                style=styles[i]
                println(io,join(("FUNCTION",h.model,h.range...,h.npx,length(h.parameters),label,
                    style.color,style.line_width,style.line_style,h.parameters...),'\t'))
                continue
            end
            if h isa ErrorGraph
                h=with_layer_style(h,label,styles[i])
                println(io,"LAYER\t",h isa BentErrorSeries ? "BENT" : "MULTI","\t",label)
                write_error_graph(io,h)
                continue
            end
            if h isa GraphSeries
                h=GraphSeries(h.x,h.y,h.exlow,h.exhigh,h.eylow,h.eyhigh,h.kind,label,h.draw,h.errors,styles[i])
                validate_graph(h)
                println(io,"LAYER\tGRAPH\t",label)
                buf=IOBuffer(); write_graph_series(buf,[h])
                for line in split(chomp(String(take!(buf))),'\n')[1:end-1] println(io,line) end
                continue
            end
            draw=get(opts,:draw,nd ? spec.draw : :hist)
            option = if draw isa AbstractString
                root_draw_option(draw)
            elseif nd
                draw==:heatmap ? "COLZ" : draw==:contour ? "CONT3" : draw==:surface ? "SURF2Z" :
                    throw(ArgumentError("unsupported ND draw preset"))
            else
                get(HIST_DRAW_OPTIONS,draw,nothing)
            end
            option === nothing && throw(ArgumentError("unsupported histogram draw preset"))
            !nd && any(s->occursin(s,uppercase(option)),("HBAR","PIE","FUNC","AXIS")) &&
                throw(ArgumentError("use standalone render_root for HBAR/PIE/FUNC/AXIS options"))
            tag=nd ? "HIST$dimensions" : "HIST1"
            println(io,join(("LAYER",tag,label,option),'\t'))
            panel_hist_style(io,styles[i])
            if h isa Hist1D
                println(io,"HIST\t",length(h.sumw))
                for (j,e) in enumerate(h.edges) println(io,"EDGE\t",j,'\t',e) end
                for j in eachindex(h.sumw) println(io,"BIN\t",j,'\t',h.sumw[j],'\t',h.sumw2[j]) end
                write_hist_stats(io,h)
            else
                values=h isa Union{Hist2D,Hist3D} ? h.sumw : h.values
                errors=h isa Union{Hist2D,Hist3D}
                println(io,join(("GRID",size(values)...,Int(errors)),'\t'))
                for (tag,edges) in (("XEDGE",h.xedges),("YEDGE",h.yedges))
                    for (j,e) in enumerate(edges) println(io,tag,'\t',j,'\t',e) end
                end
                if dimensions==3
                    for (j,e) in enumerate(h.zedges) println(io,"ZEDGE\t",j,'\t',e) end
                end
                for index in CartesianIndices(values)
                    v=errors ? h.sumw2[index] : 0.
                    println(io,join(("CELL",Tuple(index)...,values[index],v),'\t'))
                end
                errors && write_hist_stats(io,h)
            end
        end
        println(io,"END\t",length(layers))
    end
    return path
end

function write_panel_request(path,panel::PlotPanel)
    data,spec=panel.data,panel.spec
    if data isa Union{HistStack,Profile1D,Profile2D,Profile3D,PolygonHist2D}
        return write_advanced_plot(path,data,spec;panel.options...)
    elseif data isa ModelCurve
        isempty(panel.options) || throw(ArgumentError("put ModelCurve panel styling in GraphSpec"))
        return write_function_spec(path,data,spec)
    elseif data isa Union{ModelSurface,ModelVolume}
        isempty(panel.options) || throw(ArgumentError("put TF2/TF3 panel styling in Graph2DSpec/Graph3DSpec"))
        return write_analytic_spec(path,data,spec)
    elseif data isa Union{Graph3D,Graph3DErrors}
        isempty(panel.options) || throw(ArgumentError("put XYZ panel styling in Graph3DSpec"))
        validate_cloud(data); validate_nd_spec(spec)
        write_xyz_spec(path,data,spec)
        return path
    elseif data isa SpecialGraph
        isempty(panel.options) || throw(ArgumentError("put special graph styling in GraphSpec or its data style"))
        return write_special_plot(path,data,spec)
    end
    return write_basic_panel(path,panel)
end

function write_figure_grid(io,figure::PlotFigure,payloads,decorations,counter)
    figure.widths === nothing && figure.heights === nothing ||
        println(io,join(("WEIGHTS",something(figure.widths,ones(figure.layout[2]))...,
                         something(figure.heights,ones(figure.layout[1]))...),'\t'))
    for (i,cell) in enumerate(figure.panels)
        if cell isa PlotFigure
            println(io,join(("SUBGRID",i,cell.layout...,length(cell.panels)),'\t'))
            write_figure_grid(io,cell,payloads,decorations,counter)
            println(io,"ENDGRID\t",i)
        else
            k=(counter[]+=1)
            println(io,"PANEL\t",i,'\t',length(payloads[k]),'\t',length(decorations[k]))
            foreach(row->println(io,row),payloads[k])
            foreach(row->println(io,row),decorations[k])
        end
    end
end

function render_root(figure::PlotFigure,outstem::AbstractString;format=:png,root_config=nothing)
    paths=plot_output_paths(outstem,format)
    figure=figure_style(styled_figure(figure))
    leaves=figure_leaves(figure)
    mktempdir() do dir
        payloads=Vector{String}[]
        decorations=Vector{String}[]
        for (i,panel) in enumerate(leaves)
            file=write_panel_request(joinpath(dir,"panel-$i.tsv"),panel)
            push!(payloads,readlines(file)[2:end])
            length(panel.primitives)<=1000 || throw(ArgumentError("too many panel primitives"))
            push!(decorations,[primitive_record(p,panel.spec) for p in panel.primitives])
        end
        # ROOT palettes and these error/PS controls are process-global. Reject
        # incompatible panels rather than silently letting the last pad win.
        palettes=String[]
        global_controls=Tuple[]
        for (panel,rows) in zip(leaves,payloads)
            f=split(first(rows),'\t')
            f[1]=="PLOT" && parse(Int,f[4])>1 && push!(palettes,f[10])
            if hasproperty(panel.spec,:palette) push!(palettes,string(panel.spec.palette)) end
            for row in rows
                values=split(row,'\t')
                if first(values)=="CANVAS"
                    push!(global_controls,Tuple(v=="-" ? nothing : parse(Float64,v) for v in values[9:10]))
                end
            end
        end
        length(unique(palettes))<=1 || throw(ArgumentError("all pads must use the same ROOT palette"))
        # TGaxis::SetExponentOffset is static: one offset per axis direction for the whole canvas.
        for direction in ("x","y")
            offsets=Set{Tuple{String,String}}()
            for rows in payloads, row in rows
                v=split(row,'\t')
                v[1]=="AXIS" && v[2]==direction && v[end-1]!="-" && push!(offsets,(v[end-1],v[end]))
            end
            length(offsets)<=1 || throw(ArgumentError("pads have conflicting $direction exponent_offset values (TGaxis::SetExponentOffset is global)"))
        end
        for position in 1:2
            values=unique(filter(!isnothing,[c[position] for c in global_controls]))
            length(values)<=1 || throw(ArgumentError("pads have conflicting global ROOT error/PS settings"))
        end
        input=joinpath(dir,"figure.tsv")
        open(input,"w") do io
            println(io,"ROOTPlotsGraph\t1")
            println(io,join(("DRAWLAYOUT",figure.size...,figure.layout...,length(figure.panels)),'\t'))
            write_figure_grid(io,figure,payloads,decorations,Ref(0))
            println(io,"END\t",length(figure.panels))
        end
        config=find_root_config(root_config)
        config===nothing && throw(ArgumentError("local native ROOT not found"))
        exe=graphfit_executable(config)
        mkpath(dirname(outstem))
        stem=joinpath(dir,"rendered")
        run(Cmd([exe,input,stem]))
        move_plot_outputs(stem,paths)
    end
    return paths
end
render_root(panel::PlotPanel,stem::AbstractString;kwargs...) =
    render_root(PlotFigure([panel];size=panel.spec.size),stem;kwargs...)
"""
    render_root(h::Hist2D, spec::Union{Graph2DSpec,Graph3DSpec}, stem; format=:png, kwargs...)
    render_root(h::Hist3D, spec::Graph3DSpec, stem; format=:png, kwargs...)

Draw a native TH2D/TH3D (one `PlotPanel`); `kwargs` are panel keywords, e.g.
`draw="COLZ"`, `label`, `stats=true` (ROOT's Entries, Mean x/y, Std Dev x/y),
`stats_position`, `stats_name`.
"""
render_root(h::Hist2D,spec::Union{Graph2DSpec,Graph3DSpec},stem::AbstractString;
            format=:png,root_config=nothing,kwargs...) =
    render_root(PlotPanel(h,spec;kwargs...),stem;format,root_config)
render_root(h::Hist3D,spec::Graph3DSpec,stem::AbstractString;
            format=:png,root_config=nothing,kwargs...) =
    render_root(PlotPanel(h,spec;kwargs...),stem;format,root_config)
