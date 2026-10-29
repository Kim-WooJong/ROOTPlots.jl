"""Optional ROOT axis attributes. Unset fields preserve ROOT's existing defaults."""
Base.@kwdef struct AxisStyle
    title_font::Union{Nothing,Int} = nothing
    title_size::Union{Nothing,Float64} = nothing
    title_offset::Union{Nothing,Float64} = nothing
    label_font::Union{Nothing,Int} = nothing
    label_size::Union{Nothing,Float64} = nothing
    label_offset::Union{Nothing,Float64} = nothing
    tick_length::Union{Nothing,Float64} = nothing
    divisions::Union{Nothing,Int} = nothing
    max_digits::Union{Nothing,Int} = nothing
    center_title::Union{Nothing,Bool} = nothing
end

"""Optional canvas/pad and error-cap attributes for one-pad batch plots."""
Base.@kwdef struct CanvasStyle
    left_margin::Union{Nothing,Float64} = nothing
    right_margin::Union{Nothing,Float64} = nothing
    bottom_margin::Union{Nothing,Float64} = nothing
    top_margin::Union{Nothing,Float64} = nothing
    ticks_x::Union{Nothing,Bool} = nothing
    ticks_y::Union{Nothing,Bool} = nothing
    frame_line_width::Union{Nothing,Int} = nothing
    end_error_size::Union{Nothing,Float64} = nothing
    line_scale_ps::Union{Nothing,Float64} = nothing
end

"""Reference line; coordinates are physical data values, or NDC when ndc=true."""
Base.@kwdef struct PlotLine
    x1::Float64
    y1::Float64
    x2::Float64
    y2::Float64
    color::String = "#000000"
    width::Int = 2
    style::Int = 1
    ndc::Bool = false
    label::String = ""
end
PlotLine(x1,y1,x2,y2; kwargs...) = PlotLine(;x1,y1,x2,y2,kwargs...)

"""ROOT TLatex annotation; NDC coordinates by default. Uses ROOT #sigma syntax."""
Base.@kwdef struct PlotText
    x::Float64
    y::Float64
    text::String
    ndc::Bool = true
    font::Int = 42
    size::Float64 = 0.04
    color::String = "#000000"
    align::Int = 11
    angle::Float64 = 0.0
end
PlotText(x,y,text; kwargs...) = PlotText(;x,y,text=String(text),kwargs...)

"""Legend reference; :series indexes plotted data layers, :line indexes PlotLine."""
Base.@kwdef struct LegendEntry
    target::Symbol
    index::Int
    label::String
    option::String = "lep"
end
LegendEntry(target,index,label; option="lep") =
    LegendEntry(;target,index,label=String(label),option=String(option))

"""NDC legend placement. entries=nothing uses labels; [] requests an empty legend."""
Base.@kwdef struct LegendSpec
    bounds::Tuple{Float64,Float64,Float64,Float64} = (0.65,0.70,0.90,0.90)
    font::Int = 42
    text_size::Float64 = 0.04
    border_size::Int = 0
    background::Union{Nothing,String} = nothing
    columns::Int = 1
    entries::Union{Nothing,Vector{LegendEntry}} = nothing
end

valid_font(font) = 10 <= font <= 159 && mod(font,10) in (0,1,2,3)
function scene_text(value)
    occursin('\0',value) && throw(ArgumentError("NUL is not permitted in plot text"))
    clean_field(value)
    return clean_field(root_latex(value))
end
function scene_color(value)
    occursin(r"^#[0-9a-fA-F]{6}$",value) || throw(ArgumentError("color must be #RRGGBB"))
end
function validate_axis(a::AxisStyle)
    for font in (a.title_font,a.label_font)
        font === nothing || valid_font(font) || throw(ArgumentError("invalid ROOT font"))
    end
    for value in (a.title_size,a.label_size,a.title_offset)
        value === nothing || (isfinite(value) && value > 0) ||
            throw(ArgumentError("axis sizes/offset must be positive and finite"))
    end
    a.label_offset === nothing || isfinite(a.label_offset) || throw(ArgumentError("invalid label offset"))
    a.tick_length === nothing || (isfinite(a.tick_length) && 0 <= a.tick_length <= 1) ||
        throw(ArgumentError("tick length must be in 0:1"))
    a.divisions === nothing || 1 <= a.divisions <= 999999 || throw(ArgumentError("invalid axis divisions"))
    a.max_digits === nothing || 1 <= a.max_digits <= 10 || throw(ArgumentError("max_digits must be in 1:10"))
end
has_scene(spec) = spec.xaxis != AxisStyle() || spec.yaxis != AxisStyle() ||
    spec.canvas != CanvasStyle() || !isempty(spec.lines) || !isempty(spec.texts) || spec.legend_spec !== nothing

function validate_scene(spec; nseries)
    validate_axis(spec.xaxis); validate_axis(spec.yaxis)
    c=spec.canvas
    for value in (c.left_margin,c.right_margin,c.bottom_margin,c.top_margin)
        value === nothing || (isfinite(value) && 0 <= value < 1) || throw(ArgumentError("invalid canvas margin"))
    end
    # ROOT's default opposite margin is 0.1 when only one side is supplied.
    something(c.left_margin,0.1)+something(c.right_margin,0.1) < 1 &&
        something(c.bottom_margin,0.1)+something(c.top_margin,0.1) < 1 ||
        throw(ArgumentError("opposite margins must sum to less than one"))
    c.frame_line_width === nothing || 1 <= c.frame_line_width <= 100 || throw(ArgumentError("invalid frame line width"))
    c.end_error_size === nothing || (isfinite(c.end_error_size) && c.end_error_size >= 0) || throw(ArgumentError("invalid error cap size"))
    c.line_scale_ps === nothing || (isfinite(c.line_scale_ps) && c.line_scale_ps > 0) || throw(ArgumentError("invalid PS line scale"))
    length(spec.lines) <= 1000 && length(spec.texts) <= 1000 || throw(ArgumentError("too many annotations"))
    function coordinate(x,y,ndc)
        isfinite(x) && isfinite(y) || throw(ArgumentError("annotation coordinates must be finite"))
        if ndc
            0 <= x <= 1 && 0 <= y <= 1 || throw(ArgumentError("NDC coordinates must be in 0:1"))
        else
            (!spec.logx || x > 0) && (!spec.logy || y > 0) ||
                throw(ArgumentError("data annotations on log axes must be positive"))
        end
    end
    for line in spec.lines
        coordinate(line.x1,line.y1,line.ndc); coordinate(line.x2,line.y2,line.ndc)
        scene_color(line.color); scene_text(line.label)
        1 <= line.width <= 100 && 1 <= line.style <= 10 || throw(ArgumentError("invalid line style"))
    end
    for t in spec.texts
        coordinate(t.x,t.y,t.ndc); scene_color(t.color); scene_text(t.text)
        valid_font(t.font) && isfinite(t.size) && t.size > 0 &&
            div(t.align,10) in 1:3 && mod(t.align,10) in 1:3 && isfinite(t.angle) ||
            throw(ArgumentError("invalid text style"))
    end
    l=spec.legend_spec
    if l !== nothing
        x1,y1,x2,y2=l.bounds
        all(isfinite,l.bounds) && 0 <= x1 < x2 <= 1 && 0 <= y1 < y2 <= 1 || throw(ArgumentError("invalid legend bounds"))
        valid_font(l.font) && isfinite(l.text_size) && l.text_size > 0 &&
            0 <= l.border_size <= 100 && 1 <= l.columns <= 16 || throw(ArgumentError("invalid legend style"))
        l.background === nothing || scene_color(l.background)
        if l.entries !== nothing
            length(l.entries) <= 1000 || throw(ArgumentError("too many legend entries"))
            for e in l.entries
                count = e.target == :series ? nseries : e.target == :line ? length(spec.lines) : 0
                1 <= e.index <= count || throw(ArgumentError("invalid legend reference"))
                e.option in ("l","p","pe","ep","lp","lep","f") || throw(ArgumentError("invalid legend entry option"))
                scene_text(e.label)
            end
        end
    end
    return spec
end

scene_value(x) = x === nothing ? "-" : x isa Bool ? string(Int(x)) : string(x)
function write_scene(io,spec)
    has_scene(spec) || return
    println(io,"SCENE\t1")
    c=spec.canvas
    println(io,join(("CANVAS",(scene_value(getfield(c,f)) for f in fieldnames(CanvasStyle))...), '\t'))
    for (name,a) in (("x",spec.xaxis),("y",spec.yaxis))
        println(io,join(("AXIS",name,(scene_value(getfield(a,f)) for f in fieldnames(AxisStyle))...), '\t'))
    end
    for (i,l) in enumerate(spec.lines)
        println(io,join(("LINE",i,Int(l.ndc),l.x1,l.y1,l.x2,l.y2,l.color,l.width,l.style,scene_text(l.label)), '\t'))
    end
    for (i,t) in enumerate(spec.texts)
        println(io,join(("TEXT",i,Int(t.ndc),t.x,t.y,t.font,t.size,t.color,t.align,t.angle,scene_text(t.text)), '\t'))
    end
    l=spec.legend_spec
    if l !== nothing
        println(io,join(("LEGEND",l.bounds...,l.font,l.text_size,l.border_size,
                         something(l.background,"-"),l.columns,Int(l.entries !== nothing)), '\t'))
        if l.entries !== nothing
            for e in l.entries
                println(io,join(("LENTRY",e.target,e.index,e.option,scene_text(e.label)), '\t'))
            end
        end
    end
    println(io,"ENDSCENE")
end
