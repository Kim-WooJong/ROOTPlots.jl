abstract type PlotPrimitive end

"""Native ROOT arrow in physical pad coordinates."""
Base.@kwdef struct PlotArrow <: PlotPrimitive
    x1::Float64
    y1::Float64
    x2::Float64
    y2::Float64
    color::String="#000000"
    width::Int=2
    size::Float64=0.02
    option::String="|>"
end
PlotArrow(x1,y1,x2,y2;kwargs...)=PlotArrow(;x1,y1,x2,y2,kwargs...)

"""Native ROOT filled or outline rectangle in physical pad coordinates."""
Base.@kwdef struct PlotBox <: PlotPrimitive
    x1::Float64
    y1::Float64
    x2::Float64
    y2::Float64
    color::String="#000000"
    width::Int=1
    fill_color::String="#4477aa"
    fill_style::Int=0
    fill_alpha::Float64=0.3
end
PlotBox(x1,y1,x2,y2;kwargs...)=PlotBox(;x1,y1,x2,y2,kwargs...)

"""Native ROOT ellipse/arc in physical pad coordinates (angles in degrees)."""
Base.@kwdef struct PlotEllipse <: PlotPrimitive
    x::Float64
    y::Float64
    rx::Float64
    ry::Float64
    angle1::Float64=0.
    angle2::Float64=360.
    rotation::Float64=0.
    color::String="#000000"
    width::Int=1
    fill_color::String="#4477aa"
    fill_style::Int=0
    fill_alpha::Float64=0.3
end
PlotEllipse(x,y,rx,ry;kwargs...)=PlotEllipse(;x,y,rx,ry,kwargs...)

"""Native TGaxis with a linear numeric mapping; does not rescale plotted data."""
Base.@kwdef struct SecondaryAxis <: PlotPrimitive
    x1::Float64
    y1::Float64
    x2::Float64
    y2::Float64
    low::Float64
    high::Float64
    divisions::Int=510
    title::String=""
    color::String="#000000"
    option::String="+L"
end
SecondaryAxis(x1,y1,x2,y2,low,high;kwargs...)=SecondaryAxis(;x1,y1,x2,y2,low,high,kwargs...)

function primitive_record(p::PlotPrimitive,spec)
    fields=[getfield(p,k) for k in fieldnames(typeof(p))]
    all(v->!(v isa Real) || isfinite(v),fields) || throw(ArgumentError("primitive coordinates/styles must be finite"))
    scene_color(p.color)
    if p isa Union{PlotArrow,PlotBox,SecondaryAxis}
        points=((p.x1,p.y1),(p.x2,p.y2))
        (p.x1!=p.x2 || p.y1!=p.y2) || throw(ArgumentError("primitive endpoints must differ"))
    else
        p.rx>0 && p.ry>0 && p.angle1<p.angle2 || throw(ArgumentError("ellipse needs positive radii and increasing angles"))
        points=((p.x-p.rx,p.y-p.ry),(p.x+p.rx,p.y+p.ry))
    end
    for (x,y) in points
        (!spec.logx || x>0) && (!spec.logy || y>0) || throw(ArgumentError("data-coordinate primitives need positive log coordinates"))
    end
    if p isa PlotArrow
        p.width>0 && 0<p.size<=1 && p.option in (">","<","|>","<|","<>","<|>","-|") ||
            throw(ArgumentError("invalid arrow style"))
        tag="ARROW"
    elseif p isa Union{PlotBox,PlotEllipse}
        p.width>0 && 0<=p.fill_style<=9999 && 0<=p.fill_alpha<=1 || throw(ArgumentError("invalid primitive fill style"))
        scene_color(p.fill_color)
        p isa PlotBox && !(p.x1<p.x2 && p.y1<p.y2) && throw(ArgumentError("box bounds must increase"))
        tag=p isa PlotBox ? "BOX" : "ELLIPSE"
    else
        p.low<p.high && 1<=p.divisions<=999999 || throw(ArgumentError("invalid secondary axis range"))
        all(c->c in "+-LURSMWN=><",p.option) || throw(ArgumentError("unsupported secondary axis option"))
        scene_text(p.title)
        tag="AXIS"
    end
    return join(("PRIMITIVE",tag,(v isa String ? scene_text(v) : v for v in fields)...),'\t')
end
