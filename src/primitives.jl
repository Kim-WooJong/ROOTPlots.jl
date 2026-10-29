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
    p isa PlotPolygon && return polygon_record(p,spec)
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

"""
    simple_polygon(x, y) -> (x, y)

Validate one simple polygon (convex or non-convex, without holes) and return
owned vertex vectors. A repeated closing vertex equal to the first one is
removed. Requires 3:256 distinct finite vertices, no repeated neighbours, a
non-zero area, and no self-intersection (no two non-adjacent edges touch or
cross, and adjacent edges do not fold back onto each other).
"""
function simple_polygon(x::AbstractVector, y::AbstractVector)
    xx, yy = Float64.(x), Float64.(y)
    length(xx) == length(yy) || throw(ArgumentError("polygon x/y vertex arrays must have equal length"))
    if length(xx) > 3 && xx[1] == xx[end] && yy[1] == yy[end]
        pop!(xx); pop!(yy)
    end
    n = length(xx)
    3 <= n <= 256 && all(isfinite, xx) && all(isfinite, yy) ||
        throw(ArgumentError("polygon needs 3:256 finite xy vertices"))
    for i in 1:n
        j = mod1(i + 1, n)
        (xx[i] != xx[j] || yy[i] != yy[j]) || throw(ArgumentError("duplicate neighbouring polygon vertices"))
    end
    orient(ax, ay, bx, by, cx, cy) = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    area = sum(xx[i] * yy[mod1(i+1, n)] - xx[mod1(i+1, n)] * yy[i] for i in 1:n)
    isfinite(area) || throw(ArgumentError("polygon geometry overflows"))
    area != 0 || throw(ArgumentError("polygon has zero area"))
    onsegment(ax, ay, bx, by, px, py) =
        min(ax, bx) <= px <= max(ax, bx) && min(ay, by) <= py <= max(ay, by)
    function intersects(a, b, c, d)
        o1 = orient(xx[a], yy[a], xx[b], yy[b], xx[c], yy[c])
        o2 = orient(xx[a], yy[a], xx[b], yy[b], xx[d], yy[d])
        o3 = orient(xx[c], yy[c], xx[d], yy[d], xx[a], yy[a])
        o4 = orient(xx[c], yy[c], xx[d], yy[d], xx[b], yy[b])
        all(isfinite, (o1, o2, o3, o4)) || throw(ArgumentError("polygon geometry overflows"))
        (sign(o1) * sign(o2) < 0 && sign(o3) * sign(o4) < 0) && return true
        o1 == 0 && onsegment(xx[a], yy[a], xx[b], yy[b], xx[c], yy[c]) && return true
        o2 == 0 && onsegment(xx[a], yy[a], xx[b], yy[b], xx[d], yy[d]) && return true
        o3 == 0 && onsegment(xx[c], yy[c], xx[d], yy[d], xx[a], yy[a]) && return true
        o4 == 0 && onsegment(xx[c], yy[c], xx[d], yy[d], xx[b], yy[b]) && return true
        return false
    end
    for i in 1:n
        a, b, c = i, mod1(i + 1, n), mod1(i + 2, n)
        # Adjacent edges a-b and b-c must not overlap (a spike folding back).
        if orient(xx[a], yy[a], xx[b], yy[b], xx[c], yy[c]) == 0 &&
           (xx[b] - xx[a]) * (xx[c] - xx[b]) + (yy[b] - yy[a]) * (yy[c] - yy[b]) < 0
            throw(ArgumentError("polygon is self-intersecting (an edge folds back at vertex $b)"))
        end
        for j in i+2:n
            (i == 1 && j == n) && continue  # edges n-1 and 1 share vertex 1
            intersects(i, mod1(i + 1, n), j, mod1(j + 1, n)) &&
                throw(ArgumentError("polygon is self-intersecting (edges $i and $j touch or cross)"))
        end
    end
    return (xx, yy)
end

"""
    PlotPolygon(x, y; color="#000000", width=1, fill_color="#4477aa", fill_style=1001, fill_alpha=0.5)

Filled native ROOT polygon (`TPolyLine` with fill and outline) in physical pad
coordinates. Any **simple** polygon is accepted, convex or not (e.g. an L shape
or a star); see [`simple_polygon`](@ref) for the checks. `fill_style=0` draws the
outline only. Holes are not supported.
"""
struct PlotPolygon <: PlotPrimitive
    x::Vector{Float64}
    y::Vector{Float64}
    color::String
    width::Int
    fill_color::String
    fill_style::Int
    fill_alpha::Float64
    function PlotPolygon(x::AbstractVector, y::AbstractVector; color="#000000", width::Integer=1,
                         fill_color="#4477aa", fill_style::Integer=1001, fill_alpha::Real=0.5)
        xx, yy = simple_polygon(x, y)
        scene_color(color); scene_color(fill_color)
        1 <= width <= 100 && 0 <= fill_style <= 9999 && isfinite(fill_alpha) && 0 <= fill_alpha <= 1 ||
            throw(ArgumentError("invalid polygon line/fill style"))
        return new(xx, yy, String(color), Int(width), String(fill_color), Int(fill_style), Float64(fill_alpha))
    end
end

function polygon_record(p::PlotPolygon, spec)
    simple_polygon(p.x, p.y)
    for (x, y) in zip(p.x, p.y)
        (!spec.logx || x > 0) && (!spec.logy || y > 0) ||
            throw(ArgumentError("data-coordinate primitives need positive log coordinates"))
    end
    coords = Float64[]
    for i in eachindex(p.x) push!(coords, p.x[i], p.y[i]) end
    return join(("PRIMITIVE","POLYGON",length(p.x),p.color,p.width,p.fill_color,p.fill_style,p.fill_alpha,
                 coords...),'\t')
end
