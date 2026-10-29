"""Native ROOT TH1D drawing modes. `:steps` retains the original TGraph view."""
const HIST_DRAW_OPTIONS = Dict(
    :hist => "HIST", :e => "E", :e0 => "E0", :e1 => "E1",
    :e2 => "E2", :e3 => "E3", :e4 => "E4", :e5 => "E5", :e6 => "E6",
    :hist_e1 => "HIST E1", :bar => "HIST BAR",
    :points => "HIST P", :points_zero => "HIST P0",
    :line => "HIST L", :curve => "HIST C", :text => "HIST TEXT",
)
const HIST_DRAW_NAMES = Dict(value => key for (key, value) in HIST_DRAW_OPTIONS)

function hist_draw_mode(draw)
    if draw isa Symbol
        mode = draw
    elseif draw isa AbstractString
        option = uppercase(join(split(strip(draw)), " "))
        aliases = Dict("BAR" => :bar, "B" => :bar, "P" => :points,
                       "P0" => :points_zero, "L" => :line, "C" => :curve,
                       "TEXT" => :text)
        mode = option == "STEPS" ? :steps :
               get(aliases, option, get(HIST_DRAW_NAMES, option, nothing))
    else
        throw(ArgumentError("histogram draw must be a supported Symbol or ROOT option string"))
    end
    mode == :steps || haskey(HIST_DRAW_OPTIONS, mode) ||
        throw(ArgumentError("unsupported 1-D histogram draw option"))
    return mode
end

hist_error_mode(mode) = mode in (:e, :e0, :e1, :e2, :e3, :e4, :e5, :e6, :hist_e1)

"""Keep ROOT's option language opaque, but bound its C parser and TSV transport."""
function root_draw_option(option::AbstractString; allow_same=false)
    value = String(option)
    ncodeunits(value) <= 127 || throw(ArgumentError("ROOT draw option exceeds 127 bytes"))
    all(b -> 0x20 <= b <= 0x7e, codeunits(value)) ||
        throw(ArgumentError("ROOT draw option must contain printable ASCII only"))
    any(c -> c in ('"','\\',';'), value) &&
        throw(ArgumentError("quotes, backslashes, and semicolons are not valid ROOT draw options"))
    depth=0
    for c in value
        c == '(' && (depth += 1)
        c == ')' && (depth -= 1)
        depth >= 0 || throw(ArgumentError("unbalanced ROOT draw option parentheses"))
    end
    depth == 0 || throw(ArgumentError("unbalanced ROOT draw option parentheses"))
    occursin(r"GL(?:LEGO|SURF|COL|BOX|ISO)",uppercase(value)) &&
        throw(ArgumentError("OpenGL histogram options cannot be exported by the batch PNG/PDF renderer"))
    !allow_same && occursin("SAME", uppercase(value)) &&
        throw(ArgumentError("SAME needs an overlay canvas, not a standalone histogram"))
    return value
end

"""
    render_root(hists::AbstractVector{Hist1D}, spec, outstem; draw, labels, styles,
                stats=false, stats_position=:auto, stats_name=nothing, stats_headroom=false,
                format=:png)

Draw several TH1D objects on one ROOT canvas; later options must request SAME.

With `stats`, the first layer gets a ROOT statistics box (see
[`stats_option`](@ref) for the accepted values). As in ROOT, later layers get
their own box only when their option uses `SAMES` instead of `SAME`; the boxes
are stacked below the first one and their text uses the layer color. Each box
header is the layer label (the first falls back to `stats_name`).
`stats_headroom=true` raises the automatic y maximum so the data stay below the
boxes (see `render_root(::Hist1D, ::GraphSpec, stem)`).
"""
function render_root(hists::AbstractVector{Hist1D}, spec::GraphSpec, outstem::AbstractString;
                     draw, labels=fill("",length(hists)),
                     styles=fill(GraphStyle(),length(hists)), format=:png,
                     stats=false, stats_position=:auto, stats_name=nothing, stats_headroom=false,
                     root_config=nothing)
    paths = plot_output_paths(outstem, format)
    request = stats_request(; stats, stats_position, stats_headroom, ylimits=spec.ylimits,
                            name=stats_header(stats_name, isempty(labels) ? "" : first(labels)))
    2 <= length(hists) <= 64 || throw(ArgumentError("histogram overlay needs 2:64 layers"))
    length(draw) == length(labels) == length(styles) == length(hists) ||
        throw(ArgumentError("overlay options, labels, and styles must match histogram count"))
    options = String[]
    for (i, option) in enumerate(draw)
        option isa AbstractString || throw(ArgumentError("overlay draw options must be ROOT strings"))
        value = root_draw_option(option; allow_same=i>1)
        i > 1 && !occursin("SAME",uppercase(value)) &&
            throw(ArgumentError("each overlay layer after the first needs SAME"))
        push!(options,value)
    end
    limits = raw_hist_draw_limits(first(hists),spec,first(styles),first(labels); nseries=length(hists))
    for i in eachindex(hists)
        validate_stats_hist(hists[i])
        length(hists[i].sumw) <= 1_000_000 || throw(ArgumentError("too many histogram bins"))
        check_graph_style(styles[i])
        clean_field(String(labels[i]))
        spec.logx && first(hists[i].edges) <= 0 &&
            throw(ArgumentError("log x needs positive histogram edges"))
    end
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        input=joinpath(dir,"overlay.tsv")
        write_hist_overlay_spec(input,hists,spec,limits,options,labels,styles,request)
        stem=joinpath(dir,"rendered")
        run(Cmd([exe,input,stem]))
        move_plot_outputs(stem,paths)
    end
    return paths
end

function write_hist_overlay_spec(path,hists,spec,limits,options,labels,styles,request::StatsRequest)
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        println(io,join(("DRAWHISTOVERLAY",spec.size...,Int(spec.logx),Int(spec.logy),
                         limits[1]...,limits[2]...,String(spec.legend),
                         Int(spec.xlimits !== nothing),Int(spec.ylimits !== nothing),
                         Int(stats_active(request)),length(hists)),'\t'))
        write_stats_row(io,request)
        println(io,join(("LABELS",clean_field(root_latex(spec.title)),
                         clean_field(root_latex(spec.xlabel)),
                         clean_field(root_latex(spec.ylabel))),'\t'))
        write_scene(io,spec)
        for (k,h) in enumerate(hists)
            style=styles[k]
            row=("HLAYER",k,options[k],clean_field(root_latex(String(labels[k]))),
                 style.color,style.marker_style,style.marker_size,style.line_width,style.line_style)
            extended=style.fill_color!==nothing || style.fill_style!=1001 || style.fill_alpha!=0.35
            println(io,join(extended ? (row...,something(style.fill_color,style.color),style.fill_style,style.fill_alpha) : row,'\t'))
            println(io,"HIST\t",length(h.sumw))
            for (i,edge) in enumerate(h.edges) println(io,"EDGE\t",i,'\t',edge) end
            for i in eachindex(h.sumw)
                println(io,"BIN\t",i,'\t',h.sumw[i],'\t',h.sumw2[i])
            end
            write_hist_stats(io,h)
        end
        println(io,"END\t",length(hists))
    end
    return path
end

function raw_hist_draw_limits(h::Hist1D, spec::GraphSpec, style, label; nseries=1)
    validate_scene(spec; nseries)
    n = validate_stats_hist(h)
    n <= 1_000_000 || throw(ArgumentError("native histogram drawing supports at most 1,000,000 bins"))
    check_graph_style(style)
    spec.legend in (:topright, :topleft, :bottomright, :bottomleft, :none) ||
        throw(ArgumentError("invalid legend position"))
    all(s -> 1 <= s <= 10000, spec.size) ||
        throw(ArgumentError("canvas dimensions must be in 1:10000"))
    for field in (spec.title, spec.xlabel, spec.ylabel, String(label))
        clean_field(field)
        occursin('\0', field) && throw(ArgumentError("NUL is not permitted in labels"))
    end
    xlimits = spec.xlimits === nothing ? (first(h.edges), last(h.edges)) : spec.xlimits
    ylimits = spec.ylimits === nothing ? (0.0, 1.0) : spec.ylimits
    for (limits, logscale) in ((xlimits, spec.logx), (ylimits, spec.logy && spec.ylimits !== nothing))
        lo, hi = limits
        isfinite(lo) && isfinite(hi) && lo < hi && (!logscale || lo > 0) ||
            throw(ArgumentError("invalid raw histogram axis limits"))
    end
    spec.logx && first(h.edges) <= 0 &&
        throw(ArgumentError("log x needs positive histogram edges"))
    check_spec_bins(spec, "a 1-D histogram", length(h.sumw))
    return (xlimits, ylimits)
end

function hist_draw_limits(h::Hist1D, spec::GraphSpec, mode, style, label)
    validate_scene(spec; nseries=1)
    n = validate_stats_hist(h)
    n <= 1_000_000 || throw(ArgumentError("native histogram drawing supports at most 1,000,000 bins"))
    check_spec_bins(spec, "a 1-D histogram", n)
    check_graph_style(style)
    spec.legend in (:topright, :topleft, :bottomright, :bottomleft, :none) ||
        throw(ArgumentError("invalid legend position"))
    all(s -> 1 <= s <= 10000, spec.size) ||
        throw(ArgumentError("canvas dimensions must be in 1:10000"))
    for field in (spec.title, spec.xlabel, spec.ylabel, String(label))
        clean_field(field)
        occursin('\0', field) && throw(ArgumentError("NUL is not permitted in labels"))
    end
    xlo, xhi = spec.xlimits === nothing ? (first(h.edges), last(h.edges)) : spec.xlimits
    isfinite(xlo) && isfinite(xhi) && xlo < xhi &&
        (!spec.logx || (xlo > 0 && first(h.edges) > 0)) ||
        throw(ArgumentError("invalid histogram x limits for selected scale"))
    lower = copy(h.sumw)
    upper = copy(h.sumw)
    if hist_error_mode(mode)
        errors = sqrt.(h.sumw2)
        lower .-= errors
        upper .+= errors
    end
    all(isfinite, lower) && all(isfinite, upper) ||
        throw(ArgumentError("histogram error endpoints overflow Float64"))
    if spec.ylimits !== nothing
        ylo, yhi = spec.ylimits
        isfinite(ylo) && isfinite(yhi) && ylo < yhi && (!spec.logy || ylo > 0) ||
            throw(ArgumentError("invalid histogram y limits for selected scale"))
        return ((xlo, xhi), (ylo, yhi))
    end
    if spec.logy
        positive = filter(>(0), vcat(lower, upper, h.sumw))
        isempty(positive) && throw(ArgumentError("automatic log y needs positive visible data"))
        ylo = minimum(positive) / 2
        yhi = maximum(positive) * 1.46
    else
        low = min(0.0, minimum(lower))
        high = max(0.0, maximum(upper))
        if low == high
            ylo, yhi = -0.5, 0.5
        else
            pad = 0.05 * (high - low)
            ylo, yhi = low - pad, high + pad
        end
    end
    isfinite(ylo) && isfinite(yhi) && ylo < yhi && (!spec.logy || ylo > 0) ||
        throw(ArgumentError("automatic histogram y limits are not representable"))
    return ((xlo, xhi), (ylo, yhi))
end

function write_hist_draw_spec(path, h, spec, limits, mode, label, style;
                              request=stats_request(), attached_fit=nothing, overlays=GraphSeries[],
                              fit_label="")
    open(path, "w") do io
        println(io, "ROOTPlotsGraph\t1")
        if mode isa AbstractString
            println(io, join(("DRAWHISTRAW", spec.size[1], spec.size[2], Int(spec.logx), Int(spec.logy),
                              limits[1]..., limits[2]..., String(spec.legend),
                              root_draw_option(mode), Int(spec.xlimits !== nothing),
                              Int(spec.ylimits !== nothing),Int(stats_active(request))), '\t'))
        else
            println(io, join(("DRAWHIST", spec.size[1], spec.size[2], Int(spec.logx), Int(spec.logy),
                              limits[1]..., limits[2]..., String(spec.legend), String(mode)), '\t'))
        end
        write_stats_row(io, request)
        println(io, join(("LABELS", clean_field(root_latex(spec.title)),
                          clean_field(root_latex(spec.xlabel)), clean_field(root_latex(spec.ylabel)),
                          clean_field(root_latex(String(label)))), '\t'))
        write_scene(io,spec)
        base=("HSTYLE",style.color,style.marker_style,style.marker_size,style.line_width,style.line_style)
        extended=style.fill_color!==nothing || style.fill_style!=1001 || style.fill_alpha!=0.35
        println(io,join(extended ? (base...,something(style.fill_color,style.color),style.fill_style,style.fill_alpha) : base,'\t'))
        println(io, "HIST\t", length(h.sumw))
        for (i, edge) in enumerate(h.edges)
            println(io, "EDGE\t", i, '\t', edge)
        end
        for i in eachindex(h.sumw)
            println(io, "BIN\t", i, '\t', h.sumw[i], '\t', h.sumw2[i])
        end
        write_hist_stats(io, h)
        if attached_fit !== nothing
            write_attached_fit(io, attached_fit)
            isempty(fit_label) || println(io, "FITLABEL\t", clean_field(root_latex(String(fit_label))))
        end
        if isempty(overlays)
            println(io, "END\t1")
        else
            write_graph_series(io,overlays)
        end
    end
    return path
end

function render_hist_root(h::Hist1D, spec::GraphSpec, outstem::AbstractString;
                          mode, label, style, format, stats=false, fitstats=false,
                          stats_position=:auto, stats_name=nothing, fit_label="",
                          stats_headroom=false, attached_fit=nothing, overlays=GraphSeries[], root_config)
    paths = plot_output_paths(outstem, format)
    attached_fit === nothing || validate_attached_fit(attached_fit)
    request = stats_request(; stats, fitstats, stats_position, fit=attached_fit,
                            stats_headroom, ylimits=spec.ylimits,
                            name=stats_header(stats_name, label))
    clean_field(String(fit_label))
    length(overlays) <= 64 || throw(ArgumentError("too many graph overlays"))
    foreach(validate_graph,overlays)
    for s in overlays
        # Validate graph errors/centers on log axes independently of TH1 zero bins.
        graph_limits([s],with_graph_spec(spec;legend_spec=nothing,lines=PlotLine[],texts=PlotText[]))
    end
    limits = if mode isa AbstractString
        root_draw_option(mode)
        attached_fit === nothing && occursin("FUNC",uppercase(mode)) &&
            throw(ArgumentError("FUNC needs attached_fit=FitResult"))
        raw_hist_draw_limits(h, spec, style, label; nseries=1+length(overlays))
    else
        isempty(overlays) || throw(ArgumentError("mixed histogram overlays require a ROOT option string"))
        hist_draw_limits(h, spec, mode, style, label)
    end
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        input = write_hist_draw_spec(joinpath(dir, "hist.tsv"), h, spec, limits, mode, label, style;
                                     request,attached_fit,overlays,fit_label)
        stem = joinpath(dir, "rendered")
        run(Cmd([exe, input, stem]))
        move_plot_outputs(stem, paths)
    end
    return paths
end

"""Draw a native TH1D followed by graph/curve layers in one styled batch canvas.

Legend :series indices start with the histogram (1), then overlays (2,3,…).
The statistics/fit-box keywords are those of `render_root(::Hist1D, ::GraphSpec, stem)`.
"""
function render_root(h::Hist1D, overlays::AbstractVector{GraphSeries}, spec::GraphSpec,
                     outstem::AbstractString; draw="E1", label="Histogram",
                     style=GraphStyle(), format=:png, stats=false, fitstats=false,
                     stats_position=:auto, stats_name=nothing, attached_fit=nothing,
                     fit_label="", stats_headroom=false, root_config=nothing)
    draw isa AbstractString || throw(ArgumentError("mixed overlays require a ROOT draw string"))
    return render_hist_root(h,spec,outstem; mode=String(draw),label,style,format,stats,
                            fitstats,stats_position,stats_name,attached_fit,fit_label,stats_headroom,
                            overlays,root_config)
end
