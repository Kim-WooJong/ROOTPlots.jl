const ERROR_GRAPH_DRAW_MODES = (:points,:line,:linepoints,:curve,:curvepoints)

function bend_offsets(input, n, axis)
    input === nothing && return (zeros(n), zeros(n))
    if input isa AbstractVector
        lo = Float64.(input); hi = copy(lo)
    elseif input isa BinErrors || (input isa NamedTuple && keys(input) == (:low, :high))
        lo, hi = Float64.(input.low), Float64.(input.high)
    else
        throw(ArgumentError("$axis bend needs a vector or (low=…, high=…)"))
    end
    length(lo) == n == length(hi) && all(isfinite, lo) && all(isfinite, hi) ||
        throw(ArgumentError("$axis bend offsets must be finite with one value per point"))
    return (lo, hi)
end

"""
    BentErrorSeries(x, y; xerrors=nothing, yerrors=nothing, xbend=nothing, ybend=nothing,
                    label="", draw=:points, errors=:bars, style=GraphStyle())

Native ROOT `TGraphBentErrors`: asymmetric error bars whose end points can be
shifted ("bent") away from the data point, for example to separate overlapping
error bars of several series.

- `xerrors`/`yerrors` are non-negative magnitudes, as for [`GraphSeries`](@ref):
  a vector (symmetric), `(low=…, high=…)`, or `BinErrors`.
- `ybend=(low=…, high=…)` gives the **horizontal** offsets (x units) of the
  lower/upper y-error end points: the lower bar ends at `(x + low, y - eylow)`.
  A vector applies the same offset to both ends.
- `xbend=(low=…, high=…)` gives the **vertical** offsets (y units) of the
  left/right x-error end points: the left bar ends at `(x - exlow, y + low)`.
- Offsets may be negative; they are display geometry and do not change the
  stored errors.

`draw` is one of `:points`, `:line`, `:linepoints`, `:curve`, `:curvepoints`;
`errors` accepts the [`GraphSeries`](@ref) error modes. Draw it with
`render_root(series, GraphSpec(...), stem)`, in a vector of layers, or in a
[`PlotPanel`](@ref).
"""
struct BentErrorSeries
    x::Vector{Float64}
    y::Vector{Float64}
    exlow::Vector{Float64}
    exhigh::Vector{Float64}
    eylow::Vector{Float64}
    eyhigh::Vector{Float64}
    xbend_low::Vector{Float64}
    xbend_high::Vector{Float64}
    ybend_low::Vector{Float64}
    ybend_high::Vector{Float64}
    label::String
    draw::Symbol
    errors::Symbol
    style::GraphStyle
end

function BentErrorSeries(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                         xerrors=nothing, yerrors=nothing, xbend=nothing, ybend=nothing,
                         label="", draw::Symbol=:points, errors::Symbol=:bars, style=GraphStyle())
    xx, yy = Float64.(x), Float64.(y)
    n = length(xx)
    n == length(yy) && n > 0 || throw(ArgumentError("bent-error graph needs equal, nonempty x/y arrays"))
    exl, exh, _ = graph_errors(xerrors, n, :x)
    eyl, eyh, _ = graph_errors(yerrors, n, :y)
    xbl, xbh = bend_offsets(xbend, n, :x)
    ybl, ybh = bend_offsets(ybend, n, :y)
    s = BentErrorSeries(xx, yy, exl, exh, eyl, eyh, xbl, xbh, ybl, ybh, String(label), draw, errors, style)
    return validate_error_graph(s)
end

"""
    YErrorSet(errors; label="", display=:bars, style=GraphStyle(), xscale=1.0)

One y-error set of a [`MultiErrorSeries`](@ref) (for example statistical or
systematic). `errors` is a vector (symmetric), `(low=…, high=…)`, or
`BinErrors` of non-negative magnitudes.

`display` selects how this set is painted: `:bars`, `:nocaps`, `:caps`,
`:arrows`, `:filledarrows`, `:boxes` (filled rectangles), `:boxes_outline`
(filled rectangles with an outline), `:band`, `:smoothband`, or `:none`.
The line color/width/style and the fill color/style/alpha come from `style`
(`fill_color=nothing` uses `style.color`); marker fields are not used.
Boxes span the point's x errors multiplied by `xscale` (ROOT's `s=` option).
A non-empty `label` adds a legend entry for this set.
"""
struct YErrorSet
    low::Vector{Float64}
    high::Vector{Float64}
    label::String
    display::Symbol
    style::GraphStyle
    xscale::Float64
end

const ERROR_SET_DISPLAYS = (:bars,:nocaps,:caps,:arrows,:filledarrows,:boxes,:boxes_outline,
                            :band,:smoothband,:none)

function YErrorSet(errors; label="", display::Symbol=:bars, style=GraphStyle(), xscale::Real=1.0)
    errors isa YErrorSet && throw(ArgumentError("errors are already a YErrorSet"))
    if errors isa AbstractVector
        lo = Float64.(errors); hi = copy(lo)
    elseif errors isa BinErrors || (errors isa NamedTuple && keys(errors) == (:low, :high))
        lo, hi = Float64.(errors.low), Float64.(errors.high)
    else
        throw(ArgumentError("y-error set needs a vector, (low=…, high=…), or BinErrors"))
    end
    return YErrorSet(lo, hi, String(label), display, style, Float64(xscale))
end

"""
    MultiErrorSeries(x, y, sets; xerrors=nothing, label="", draw=:points,
                     sum=:quadrature, style=GraphStyle())

Native ROOT `TGraphMultiErrors`: one set of points with several y-error sets,
e.g. statistical and systematic uncertainties, each painted with its own style.
`sets` is a vector of [`YErrorSet`](@ref) values; plain error vectors or
`(low=…, high=…)` tuples are wrapped with default settings. At most 16 sets
are accepted.

`sum` is ROOT's error-summation mode (`:first` = `kOnlyFirst`, `:quadrature` =
`kSquareSum`, `:linear` = `kAbsSum`). It is stored on the native object and
used by [`total_yerrors`](@ref); ROOT paints every set separately from the
point, so the mode does **not** change the drawing. Whether a combination rule
is appropriate (independence, correlations) is the caller's scientific decision.

`draw` is `:points`, `:line`, `:linepoints`, `:curve`, or `:curvepoints`.
`style` sets markers and connecting lines; x errors (if any) use its line style.
"""
struct MultiErrorSeries
    x::Vector{Float64}
    y::Vector{Float64}
    exlow::Vector{Float64}
    exhigh::Vector{Float64}
    sets::Vector{YErrorSet}
    label::String
    draw::Symbol
    sum::Symbol
    style::GraphStyle
end

function MultiErrorSeries(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, sets::AbstractVector;
                          xerrors=nothing, label="", draw::Symbol=:points, sum::Symbol=:quadrature,
                          style=GraphStyle())
    xx, yy = Float64.(x), Float64.(y)
    n = length(xx)
    n == length(yy) && n > 0 || throw(ArgumentError("multi-error graph needs equal, nonempty x/y arrays"))
    exl, exh, _ = graph_errors(xerrors, n, :x)
    owned = YErrorSet[]
    for set in sets
        set isa YErrorSet || (set = YErrorSet(set))
        push!(owned, YErrorSet(copy(set.low), copy(set.high), set.label, set.display, set.style, set.xscale))
    end
    s = MultiErrorSeries(xx, yy, exl, exh, owned, String(label), draw, sum, style)
    return validate_error_graph(s)
end

const ErrorGraph = Union{BentErrorSeries,MultiErrorSeries}

function validate_error_graph(s::BentErrorSeries)
    n = length(s.x)
    n > 0 && all(length(v) == n for v in (s.y,s.exlow,s.exhigh,s.eylow,s.eyhigh,
                                           s.xbend_low,s.xbend_high,s.ybend_low,s.ybend_high)) ||
        throw(ArgumentError("bent-error graph arrays changed length"))
    all(isfinite, s.x) && all(isfinite, s.y) || throw(ArgumentError("graph coordinates must be finite"))
    for e in (s.exlow,s.exhigh,s.eylow,s.eyhigh)
        all(isfinite, e) && all(>=(0), e) || throw(ArgumentError("bent-error magnitudes must be finite and non-negative"))
    end
    for b in (s.xbend_low,s.xbend_high,s.ybend_low,s.ybend_high)
        all(isfinite, b) || throw(ArgumentError("bend offsets must be finite"))
    end
    for v in (s.x .- s.exlow, s.x .+ s.exhigh, s.y .- s.eylow, s.y .+ s.eyhigh,
              s.y .+ s.xbend_low, s.y .+ s.xbend_high, s.x .+ s.ybend_low, s.x .+ s.ybend_high)
        all(isfinite, v) || throw(ArgumentError("bent-error end points overflow"))
    end
    s.draw in ERROR_GRAPH_DRAW_MODES || throw(ArgumentError("bent-error draw must be one of $(ERROR_GRAPH_DRAW_MODES)"))
    s.errors in GRAPH_ERROR_MODES || throw(ArgumentError("invalid bent-error display mode"))
    scene_text(s.label)
    check_graph_style(s.style)
    return s
end

function validate_error_graph(s::MultiErrorSeries)
    n = length(s.x)
    n > 0 && length(s.y) == length(s.exlow) == length(s.exhigh) == n ||
        throw(ArgumentError("multi-error graph arrays changed length"))
    1 <= length(s.sets) <= 16 || throw(ArgumentError("multi-error graph needs 1:16 y-error sets"))
    all(isfinite, s.x) && all(isfinite, s.y) || throw(ArgumentError("graph coordinates must be finite"))
    for e in (s.exlow, s.exhigh)
        all(isfinite, e) && all(>=(0), e) || throw(ArgumentError("x errors must be finite and non-negative"))
    end
    for (k, set) in enumerate(s.sets)
        length(set.low) == length(set.high) == n ||
            throw(ArgumentError("y-error set $k needs one low/high value per point"))
        all(isfinite, set.low) && all(isfinite, set.high) && all(>=(0), set.low) && all(>=(0), set.high) ||
            throw(ArgumentError("y-error set $k must be finite and non-negative"))
        all(isfinite, s.y .- set.low) && all(isfinite, s.y .+ set.high) ||
            throw(ArgumentError("y-error set $k end points overflow"))
        set.display in ERROR_SET_DISPLAYS ||
            throw(ArgumentError("y-error set display must be one of $(ERROR_SET_DISPLAYS)"))
        isfinite(set.xscale) && 0 < set.xscale <= 100 ||
            throw(ArgumentError("y-error set xscale must be in (0, 100]"))
        scene_text(set.label)
        check_graph_style(set.style)
    end
    all(isfinite, s.x .- s.exlow) && all(isfinite, s.x .+ s.exhigh) ||
        throw(ArgumentError("x error end points overflow"))
    s.draw in ERROR_GRAPH_DRAW_MODES || throw(ArgumentError("multi-error draw must be one of $(ERROR_GRAPH_DRAW_MODES)"))
    s.sum in (:first, :quadrature, :linear) || throw(ArgumentError("sum must be :first, :quadrature, or :linear"))
    scene_text(s.label)
    check_graph_style(s.style)
    return s
end

"""
    total_yerrors(g::MultiErrorSeries) -> (low=…, high=…)

Combine the y-error sets with the series' `sum` rule, as ROOT's
`TGraphMultiErrors::GetEYlow/GetEYhigh` do: `:first` returns the first set,
`:quadrature` the square root of the summed squares, `:linear` the plain sum.
The result can, e.g., be passed as `yerrors` of a `GraphSeries` for `fit_root`.
"""
function total_yerrors(g::MultiErrorSeries)
    validate_error_graph(g)
    g.sum == :first && return (low=copy(g.sets[1].low), high=copy(g.sets[1].high))
    if g.sum == :quadrature
        return (low=sqrt.(sum(s.low.^2 for s in g.sets)), high=sqrt.(sum(s.high.^2 for s in g.sets)))
    end
    return (low=sum(s.low for s in g.sets), high=sum(s.high for s in g.sets))
end

# A plain GraphSeries whose error envelope covers everything ROOT paints; used
# only for frame ranges and log-axis validation.
function error_graph_envelope(s::BentErrorSeries)
    s.errors == :none && return GraphSeries(s.x, s.y; draw=s.draw, errors=:none, style=s.style)
    left = min.(s.x .- s.exlow, s.x .+ s.ybend_low, s.x .+ s.ybend_high)
    right = max.(s.x .+ s.exhigh, s.x .+ s.ybend_low, s.x .+ s.ybend_high)
    bottom = min.(s.y .- s.eylow, s.y .+ s.xbend_low, s.y .+ s.xbend_high)
    top = max.(s.y .+ s.eyhigh, s.y .+ s.xbend_low, s.y .+ s.xbend_high)
    return GraphSeries(s.x, s.y; xerrors=(low=s.x .- left, high=right .- s.x),
                       yerrors=(low=s.y .- bottom, high=top .- s.y), draw=s.draw, style=s.style)
end

function error_graph_envelope(s::MultiErrorSeries)
    shown = [set for set in s.sets if set.display != :none]
    n = length(s.x)
    low = isempty(shown) ? zeros(n) : reduce((a, b) -> max.(a, b), [set.low for set in shown])
    high = isempty(shown) ? zeros(n) : reduce((a, b) -> max.(a, b), [set.high for set in shown])
    return GraphSeries(s.x, s.y; xerrors=(low=copy(s.exlow), high=copy(s.exhigh)),
                       yerrors=(low=low, high=high), draw=s.draw, style=s.style)
end

error_graph_limits(s::ErrorGraph, scene) = graph_limits([error_graph_envelope(validate_error_graph(s))], scene)

with_layer_style(s::BentErrorSeries, label, style) =
    validate_error_graph(BentErrorSeries(s.x,s.y,s.exlow,s.exhigh,s.eylow,s.eyhigh,s.xbend_low,s.xbend_high,
                                         s.ybend_low,s.ybend_high,String(label),s.draw,s.errors,style))
with_layer_style(s::MultiErrorSeries, label, style) =
    validate_error_graph(MultiErrorSeries(s.x,s.y,s.exlow,s.exhigh,s.sets,String(label),s.draw,s.sum,style))

function write_error_graph(io, s::BentErrorSeries)
    st = s.style
    println(io, join(("BENT",length(s.x),s.draw,s.errors,scene_text(s.label),st.color,st.marker_style,
                      st.marker_size,st.line_width,st.line_style,something(st.fill_color,"-"),
                      st.fill_style,st.fill_alpha),'\t'))
    for i in eachindex(s.x)
        println(io, join(("BPOINT",i,s.x[i],s.y[i],s.exlow[i],s.exhigh[i],s.eylow[i],s.eyhigh[i],
                          s.xbend_low[i],s.xbend_high[i],s.ybend_low[i],s.ybend_high[i]),'\t'))
    end
end

function write_error_graph(io, s::MultiErrorSeries)
    st = s.style
    xerrors = any(>(0), s.exlow) || any(>(0), s.exhigh)
    println(io, join(("MULTI",length(s.x),length(s.sets),s.draw,s.sum,scene_text(s.label),st.color,
                      st.marker_style,st.marker_size,st.line_width,st.line_style,Int(xerrors)),'\t'))
    for (k, set) in enumerate(s.sets)
        e = set.style
        println(io, join(("ESET",k,set.display,scene_text(set.label),e.color,e.line_width,e.line_style,
                          something(e.fill_color,e.color),e.fill_style,e.fill_alpha,set.xscale),'\t'))
    end
    for i in eachindex(s.x)
        values = Float64[]
        for set in s.sets push!(values, set.low[i], set.high[i]) end
        println(io, join(("MPOINT",i,s.x[i],s.y[i],s.exlow[i],s.exhigh[i],values...),'\t'))
    end
end

render_root(s::ErrorGraph, spec::GraphSpec, outstem::AbstractString; format=:png, root_config=nothing) =
    render_root([s], spec, outstem; format, root_config)

"""
    render_root(layers::AbstractVector, spec::GraphSpec, outstem; format=:png, root_config=nothing)

Draw mixed XY layers (`GraphSeries`, `BentErrorSeries`, `MultiErrorSeries`,
`ModelCurve`, `Hist1D`) on one pad. Equivalent to
`render_root(PlotPanel(layers, spec), outstem)`; vectors containing only
`GraphSeries` or only `Hist1D` keep their dedicated methods.
"""
function render_root(layers::AbstractVector, spec::GraphSpec, outstem::AbstractString;
                     format=:png, root_config=nothing)
    panel_layers(layers)  # rejects unsupported layer types before any output
    return render_root(PlotFigure([PlotPanel(collect(layers), spec)]; size=spec.size), outstem; format, root_config)
end
