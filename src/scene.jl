"""
    LabelChange(index; angle=nothing, size=nothing, align=nothing, color=nothing,
                font=nothing, text=nothing)

One ROOT `TAxis::ChangeLabel(labNum, labAngle, labSize, labAlign, labColor,
labFont, labText)` call for a numeric axis label. `index` is ROOT's label number:
1 is the first drawn label, -1 the last (negative numbers count from the end);
0 is invalid. Unset fields keep the label's attribute (ROOT's -1 / empty text).
`size=0` hides the label; `color` is `#RRGGBB`; `align` is a ROOT text
alignment (11:33); `angle` is in degrees; `text` replaces the label (ROOT TLatex).
Labels of alphanumeric (bin-label) axes are not affected, as in ROOT.
"""
struct LabelChange
    index::Int
    angle::Union{Nothing,Float64}
    size::Union{Nothing,Float64}
    align::Union{Nothing,Int}
    color::Union{Nothing,String}
    font::Union{Nothing,Int}
    text::Union{Nothing,String}
end
LabelChange(index::Integer; angle=nothing, size=nothing, align=nothing, color=nothing,
            font=nothing, text=nothing) =
    LabelChange(Int(index), angle, size, align, color, font, text)

"""
    AxisStyle(; kwargs...)

Optional ROOT axis attributes, applied with the named ROOT setter to the
drawn axis (`TAxis`; palette axes import the histogram's z axis). Every field
defaults to `nothing`, which keeps ROOT's default (or the package's
documented layout default) unchanged.

| field | ROOT call |
|---|---|
| `title_font`, `title_size`, `title_offset`, `title_color` | `TAxis::SetTitleFont/Size/Offset/Color` |
| `label_font`, `label_size`, `label_offset`, `label_color` | `TAxis::SetLabelFont/Size/Offset/Color` |
| `axis_color` | `TAxis::SetAxisColor` |
| `tick_length` | `TAxis::SetTickLength` |
| `ticks` (`"+"`, `"-"`, `"+-"`) | `TAxis::SetTicks` (tick side; see `CanvasStyle.ticks_x` for the opposite frame edge) |
| `divisions`, `optimize_divisions` | `TAxis::SetNdivisions(n, optim)` (`optimize_divisions=false` keeps exactly `n`) |
| `max_digits` | `TAxis::SetMaxDigits` |
| `no_exponent` | `TAxis::SetNoExponent` (print 10000 instead of 10^{4}) |
| `more_log_labels` | `TAxis::SetMoreLogLabels` (label intermediate log ticks) |
| `decimals` | `TAxis::SetDecimals` (same number of decimals on every label) |
| `center_title`, `center_labels` | `TAxis::CenterTitle`, `TAxis::CenterLabels` |
| `time_display` | `TAxis::SetTimeDisplay` |
| `time_format` | `TAxis::SetTimeFormat` (strftime format; ROOT's `%F` offset syntax allowed) |
| `time_offset` | `TAxis::SetTimeOffset(t, "gmt")`: UNIX time (s) of axis value 0, in GMT |
| `exponent_offset` `(dx, dy)` | `TGaxis::SetExponentOffset(dx, dy, axis)` with axis `"x"` or `"y"` (static; one value per direction per canvas) |
| `bin_labels` | `TAxis::SetBinLabel(bin, label)` (histogram bin axes only) |
| `labels_option` (`"h"`, `"v"`, `"u"`, `"d"`) | `TAxis::LabelsOption` (orientation of bin labels; needs `bin_labels`) |
| `change_labels` | `TAxis::ChangeLabel` per [`LabelChange`](@ref) |

`bin_labels` accepts a vector of strings (bins 1, 2, …), `bin => label` pairs, or
a dictionary. Precision-3 fonts (e.g. 43) use pixel sizes, as in ROOT.
"""
struct AxisStyle
    title_font::Union{Nothing,Int}
    title_size::Union{Nothing,Float64}
    title_offset::Union{Nothing,Float64}
    label_font::Union{Nothing,Int}
    label_size::Union{Nothing,Float64}
    label_offset::Union{Nothing,Float64}
    tick_length::Union{Nothing,Float64}
    divisions::Union{Nothing,Int}
    max_digits::Union{Nothing,Int}
    center_title::Union{Nothing,Bool}
    optimize_divisions::Union{Nothing,Bool}
    no_exponent::Union{Nothing,Bool}
    more_log_labels::Union{Nothing,Bool}
    decimals::Union{Nothing,Bool}
    center_labels::Union{Nothing,Bool}
    label_color::Union{Nothing,String}
    title_color::Union{Nothing,String}
    axis_color::Union{Nothing,String}
    ticks::Union{Nothing,String}
    labels_option::Union{Nothing,String}
    time_display::Union{Nothing,Bool}
    time_format::Union{Nothing,String}
    time_offset::Union{Nothing,Float64}
    exponent_offset::Union{Nothing,Tuple{Float64,Float64}}
    bin_labels::Union{Nothing,Vector{Pair{Int,String}}}
    change_labels::Union{Nothing,Vector{LabelChange}}
end

function AxisStyle(; title_font=nothing, title_size=nothing, title_offset=nothing,
                   label_font=nothing, label_size=nothing, label_offset=nothing,
                   tick_length=nothing, divisions=nothing, max_digits=nothing,
                   center_title=nothing, optimize_divisions=nothing, no_exponent=nothing,
                   more_log_labels=nothing, decimals=nothing, center_labels=nothing,
                   label_color=nothing, title_color=nothing, axis_color=nothing,
                   ticks=nothing, labels_option=nothing, time_display=nothing,
                   time_format=nothing, time_offset=nothing, exponent_offset=nothing,
                   bin_labels=nothing, change_labels=nothing)
    return AxisStyle(title_font, title_size, title_offset, label_font, label_size, label_offset,
                     tick_length, divisions, max_digits, center_title, optimize_divisions,
                     no_exponent, more_log_labels, decimals, center_labels, label_color,
                     title_color, axis_color, ticks === nothing ? nothing : String(ticks),
                     labels_option === nothing ? nothing : String(labels_option),
                     time_display, time_format, time_offset,
                     exponent_offset === nothing ? nothing : Tuple(Float64.(exponent_offset)),
                     axis_bin_labels(bin_labels), axis_label_changes(change_labels))
end

# Vector-valued fields compare by value, so AxisStyle() == AxisStyle() and equal styles hash equally.
Base.:(==)(a::AxisStyle, b::AxisStyle) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(AxisStyle))
Base.hash(a::AxisStyle, h::UInt) = hash(Tuple(getfield(a, f) for f in fieldnames(AxisStyle)), hash(:AxisStyle, h))

axis_bin_labels(::Nothing) = nothing
function axis_bin_labels(labels)
    labels isa Union{AbstractString,Symbol,Number} &&
        throw(ArgumentError("bin_labels needs a vector of labels, bin => label pairs, or a Dict"))
    out = Pair{Int,String}[]
    items = labels isa AbstractDict ? collect(labels) : collect(labels)
    if !isempty(items) && all(x -> x isa AbstractString, items)
        for (i, label) in enumerate(items) push!(out, i => String(label)) end
        return out
    end
    for item in items
        item isa Pair && first(item) isa Integer && last(item) isa AbstractString ||
            throw(ArgumentError("bin_labels needs a vector of labels, bin => label pairs, or a Dict"))
        push!(out, Int(first(item)) => String(last(item)))
    end
    return sort!(out; by=first)
end
axis_label_changes(::Nothing) = nothing
axis_label_changes(c::LabelChange) = [c]
function axis_label_changes(changes)
    all(x -> x isa LabelChange, changes) ||
        throw(ArgumentError("change_labels needs LabelChange values"))
    return LabelChange[changes...]
end

"""
    CanvasStyle(; kwargs...)

Optional pad (`TPad`) and pad-title attributes. Every field defaults to `nothing`
(ROOT's or the package's documented default).

| field | ROOT call |
|---|---|
| `left_margin`, `right_margin`, `bottom_margin`, `top_margin` | `TPad::SetLeftMargin` … (pad fractions) |
| `ticks_x`, `ticks_y` | `TPad::SetTickx/SetTicky` (ticks on the opposite frame edge) |
| `grid_x`, `grid_y` | `TPad::SetGridx/SetGridy` |
| `frame_line_width` | `TPad::SetFrameLineWidth` |
| `end_error_size`, `line_scale_ps` | `gStyle->SetEndErrorSize`, `gStyle->SetLineScalePS` (process-global) |
| `opt_title` | `gStyle->SetOptTitle`: `true` lets ROOT draw the plot title as its `TPaveText` "title" box |
| `title_x`, `title_y`, `title_w`, `title_h` | `gStyle->SetTitleX/Y/W/H` (NDC; 0 = ROOT's automatic size) |
| `title_align` | `gStyle->SetTitleAlign` (e.g. 23 = centred at the top, ROOT's default) |
| `title_font`, `title_font_size` | `gStyle->SetTitleFont(font, "t")`, `gStyle->SetTitleFontSize` |
| `title_border_size`, `title_fill_color`, `title_style`, `title_text_color` | `gStyle->SetTitleBorderSize/SetTitleFillColor/SetTitleStyle/SetTitleTextColor` |

Without `opt_title=true` (the default) a non-empty plot title is drawn as the
package's TLatex heading at the top left, as before. With `opt_title=true` the
title string is passed to `TH1::SetTitle`, so ROOT's `"title;x title;y title"`
form also sets the axis titles. The `title_*` fields require `opt_title=true`;
they apply to this pad only (the values are set while ROOT creates the box).
"""
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
    grid_x::Union{Nothing,Bool} = nothing
    grid_y::Union{Nothing,Bool} = nothing
    opt_title::Union{Nothing,Bool} = nothing
    title_x::Union{Nothing,Float64} = nothing
    title_y::Union{Nothing,Float64} = nothing
    title_w::Union{Nothing,Float64} = nothing
    title_h::Union{Nothing,Float64} = nothing
    title_align::Union{Nothing,Int} = nothing
    title_font::Union{Nothing,Int} = nothing
    title_font_size::Union{Nothing,Float64} = nothing
    title_border_size::Union{Nothing,Int} = nothing
    title_fill_color::Union{Nothing,String} = nothing
    title_style::Union{Nothing,Int} = nothing
    title_text_color::Union{Nothing,String} = nothing
end

const CANVAS_TITLE_FIELDS = (:title_x,:title_y,:title_w,:title_h,:title_align,:title_font,
                             :title_font_size,:title_border_size,:title_fill_color,:title_style,
                             :title_text_color)

"""Overlay the explicitly set fields of `over` on `base` (layout defaults)."""
merge_canvas(base::CanvasStyle, over::CanvasStyle) =
    CanvasStyle(; (f => something(getfield(over, f), Some(getfield(base, f))) for f in fieldnames(CanvasStyle))...)

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
printable_ascii(value) = all(b -> 0x20 <= b <= 0x7e, codeunits(value))

function validate_label_change(c::LabelChange)
    c.index != 0 && abs(c.index) <= 1000 ||
        throw(ArgumentError("LabelChange index must be nonzero (1 = first label, -1 = last) with |index| <= 1000"))
    c.angle === nothing || (isfinite(c.angle) && -360 <= c.angle <= 360) ||
        throw(ArgumentError("LabelChange angle must be finite degrees in -360:360"))
    c.size === nothing || (isfinite(c.size) && c.size >= 0) ||
        throw(ArgumentError("LabelChange size must be finite and nonnegative (0 hides the label)"))
    c.align === nothing || (div(c.align,10) in 1:3 && mod(c.align,10) in 1:3) ||
        throw(ArgumentError("LabelChange align must be a ROOT text alignment 11:33"))
    c.color === nothing || scene_color(c.color)
    c.font === nothing || valid_font(c.font) || throw(ArgumentError("invalid LabelChange font"))
    if c.text !== nothing
        isempty(c.text) && throw(ArgumentError("LabelChange text must not be empty; use size=0 to hide a label"))
        scene_text(c.text)
    end
    return c
end

"""Validate an AxisStyle for axis `name` (:x, :y, :z, or a two-pad axis such as :ratio)."""
function validate_axis(a::AxisStyle; name=:x)
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
    for color in (a.label_color,a.title_color,a.axis_color)
        color === nothing || scene_color(color)
    end
    a.ticks === nothing || a.ticks in ("+","-","+-") ||
        throw(ArgumentError("ticks must be \"+\", \"-\", or \"+-\" (TAxis::SetTicks)"))
    a.labels_option === nothing || a.labels_option in ("h","v","u","d") ||
        throw(ArgumentError("labels_option must be \"h\", \"v\", \"u\", or \"d\" (TAxis::LabelsOption)"))
    a.labels_option !== nothing && a.bin_labels === nothing &&
        throw(ArgumentError("labels_option orients bin labels and needs bin_labels (TAxis::LabelsOption)"))
    if a.time_format !== nothing
        f=a.time_format
        1 <= ncodeunits(f) <= 256 && printable_ascii(f) ||
            throw(ArgumentError("time_format must be 1:256 printable ASCII characters"))
    end
    a.time_offset === nothing || isfinite(a.time_offset) ||
        throw(ArgumentError("time_offset must be a finite UNIX time in seconds"))
    if a.exponent_offset !== nothing
        name in (:x,:y,:ratio) ||
            throw(ArgumentError("exponent_offset applies to x/y axes only (TGaxis::SetExponentOffset)"))
        all(v -> isfinite(v) && abs(v) <= 1, a.exponent_offset) ||
            throw(ArgumentError("exponent_offset needs two finite NDC offsets in -1:1"))
    end
    if a.bin_labels !== nothing
        1 <= length(a.bin_labels) <= 10_000 ||
            throw(ArgumentError("bin_labels needs 1:10000 labels; use nothing for none"))
        allunique(first.(a.bin_labels)) || throw(ArgumentError("bin_labels must name each bin once"))
        for (bin,label) in a.bin_labels
            1 <= bin <= 1_000_000 || throw(ArgumentError("bin label index must be a histogram bin number >= 1"))
            isempty(label) && throw(ArgumentError("bin labels must not be empty"))
            scene_text(label)
        end
    end
    if a.change_labels !== nothing
        length(a.change_labels) <= 1000 || throw(ArgumentError("too many change_labels"))
        foreach(validate_label_change, a.change_labels)
    end
    return a
end

"""
Check `bin_labels` against the bins of the drawn histogram axis: `nbins` is the
number of bins, or `nothing` when the axis is not a histogram bin axis.
"""
function check_axis_bins(a::AxisStyle, nbins, name, what)
    a.bin_labels === nothing && return a
    nbins === nothing &&
        throw(ArgumentError("$name-axis bin_labels need a histogram bin axis, which $what does not draw " *
                            "(TAxis::SetBinLabel); use change_labels to relabel numeric axes"))
    maximum(first, a.bin_labels) <= nbins ||
        throw(ArgumentError("$name-axis bin label index exceeds the $nbins histogram bins"))
    return a
end
"""Bin-label checks for the x/y(/z) AxisStyles of a spec; `bins` holds one entry per axis."""
function check_spec_bins(spec, what, bins...)
    for ((name,a),n) in zip(scene_axes(spec), (bins..., nothing, nothing, nothing))
        check_axis_bins(a, n, name, what)
    end
    return spec
end

function validate_canvas(c::CanvasStyle)
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
    if c.opt_title !== true
        unused = [f for f in CANVAS_TITLE_FIELDS if getfield(c,f) !== nothing]
        isempty(unused) || throw(ArgumentError("$(join(unused, ", ")) need opt_title=true (ROOT title box)"))
    end
    for value in (c.title_x,c.title_y,c.title_w,c.title_h)
        value === nothing || (isfinite(value) && 0 <= value <= 1) ||
            throw(ArgumentError("title_x/title_y/title_w/title_h must be NDC values in 0:1"))
    end
    c.title_align === nothing || (div(c.title_align,10) in 1:3 && mod(c.title_align,10) in 1:3) ||
        throw(ArgumentError("title_align must be a ROOT text alignment 11:33"))
    c.title_font === nothing || valid_font(c.title_font) || throw(ArgumentError("invalid title font"))
    c.title_font_size === nothing || (isfinite(c.title_font_size) && 0 < c.title_font_size <= 1000) ||
        throw(ArgumentError("title_font_size must be positive (NDC, or pixels for precision-3 fonts)"))
    c.title_border_size === nothing || 0 <= c.title_border_size <= 100 ||
        throw(ArgumentError("title_border_size must be in 0:100"))
    c.title_style === nothing || c.title_style in (0,1001) || 3000 <= c.title_style <= 4100 ||
        throw(ArgumentError("title_style must be a ROOT fill style (0, 1001, 3000:4100)"))
    for color in (c.title_fill_color,c.title_text_color)
        color === nothing || scene_color(color)
    end
    return c
end

scene_axes(spec) = hasproperty(spec,:zaxis) ? ((:x,spec.xaxis),(:y,spec.yaxis),(:z,spec.zaxis)) :
                                              ((:x,spec.xaxis),(:y,spec.yaxis))
has_scene(spec) = spec.xaxis != AxisStyle() || spec.yaxis != AxisStyle() ||
    spec.canvas != CanvasStyle() || !isempty(spec.lines) || !isempty(spec.texts) || spec.legend_spec !== nothing ||
    has_style_or_label(spec)

function validate_scene(spec; nseries)
    validate_axis(spec.xaxis; name=:x); validate_axis(spec.yaxis; name=:y)
    validate_canvas(spec.canvas)
    validate_style_and_label(spec)
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
const AXIS_TICKS_CODE = Dict("+"=>"plus", "-"=>"minus", "+-"=>"both")
"""Scalar AxisStyle fields in protocol order ("-" = unset)."""
function axis_scalars(a::AxisStyle)
    ex = a.exponent_offset
    return (scene_value(a.title_font), scene_value(a.title_size), scene_value(a.title_offset),
            scene_value(a.label_font), scene_value(a.label_size), scene_value(a.label_offset),
            scene_value(a.tick_length), scene_value(a.divisions), scene_value(a.max_digits),
            scene_value(a.center_title), scene_value(a.optimize_divisions), scene_value(a.no_exponent),
            scene_value(a.more_log_labels), scene_value(a.decimals), scene_value(a.center_labels),
            scene_value(a.label_color), scene_value(a.title_color), scene_value(a.axis_color),
            a.ticks === nothing ? "-" : AXIS_TICKS_CODE[a.ticks], scene_value(a.labels_option),
            scene_value(a.time_display), scene_value(a.time_offset),
            ex === nothing ? "-" : string(ex[1]), ex === nothing ? "-" : string(ex[2]))
end
"""One axis record (`tag`, axis name, scalars) followed by its time-format, bin-label and label-change rows."""
function write_axis(io, tag, name, a::AxisStyle)
    println(io, join((tag, name, axis_scalars(a)...), '\t'))
    a.time_format === nothing || println(io, join(("AXTIME", name, clean_field(a.time_format)), '\t'))
    if a.bin_labels !== nothing
        for (bin,label) in a.bin_labels
            println(io, join(("AXBIN", name, bin, scene_text(label)), '\t'))
        end
    end
    if a.change_labels !== nothing
        for c in a.change_labels
            println(io, join(("AXCHANGE", name, c.index, something(c.angle,-1.0), something(c.size,-1.0),
                              something(c.align,-1), something(c.color,"-"), something(c.font,-1),
                              Int(c.text !== nothing), c.text === nothing ? "-" : scene_text(c.text)), '\t'))
        end
    end
end
function write_scene(io,spec)
    has_scene(spec) || return
    println(io,"SCENE\t1")
    spec.style === nothing || println(io,"ROOTSTYLE\t",String(spec.style))
    c=spec.canvas
    println(io,join(("CANVAS",(scene_value(getfield(c,f)) for f in fieldnames(CanvasStyle))...), '\t'))
    for (name,a) in scene_axes(spec)
        write_axis(io,"AXIS",String(name),a)
    end
    hasproperty(spec,:lines) || (write_experiment_label(io,experiment_label(spec)); println(io,"ENDSCENE"); return)
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
    write_experiment_label(io,experiment_label(spec))
    println(io,"ENDSCENE")
end

# ---------------------------------------------------------------------------
# ROOT styles (gROOT->SetStyle) and experiment labels (ROOT TLatex)

const LABEL_POSITIONS = (:inside_left, :inside_right, :outside)
const CMS_LABEL_KEYS = (:extra, :lumi, :energy, :position)
const ATLAS_LABEL_KEYS = (:text, :energy, :lumi, :position)

"""
Check a `style` value: `nothing` (ROOT's default gStyle, untouched) or a Symbol
naming a ROOT `TStyle`. The native backend checks the name against gROOT's
style list (see [`root_styles`](@ref)); names match case-insensitively.
"""
function check_style_name(style)
    style === nothing && return nothing
    style isa Symbol ||
        throw(ArgumentError("style must be nothing or a Symbol naming a ROOT TStyle, e.g. :CMS, :ATLAS, :Plain"))
    name = String(style)
    1 <= ncodeunits(name) <= 64 && occursin(r"^[A-Za-z0-9_-]+$", name) ||
        throw(ArgumentError("invalid ROOT style name $(repr(style))"))
    return style
end
style_is(style, name) = style !== nothing && lowercase(String(style)) == lowercase(name)

"""
Normalized experiment label of one pad: `kind` ("CMS" or "ATLAS"), `position`,
the text after/below the experiment name (`extra`), the text above the frame on
the right (`right`, CMS luminosity/energy), and further lines below the label.
"""
struct ExperimentLabel
    kind::String
    position::Symbol
    extra::String
    right::String
    lines::Vector{String}
end

function label_options(value, keys, name)
    value === true && return NamedTuple()
    value isa NamedTuple ||
        throw(ArgumentError("$name must be nothing, true, false, or a NamedTuple with keys $(join(keys, ", "))"))
    unknown = [k for k in Base.keys(value) if !(k in keys)]
    isempty(unknown) || throw(ArgumentError("unknown $name keys $(join(unknown, ", ")); use $(join(keys, ", "))"))
    return value
end
function label_string(o, key, default, name)
    v = get(o, key, default)
    v isa AbstractString || throw(ArgumentError("$name.$key must be a string"))
    scene_text(String(v))
    ncodeunits(v) <= 256 || throw(ArgumentError("$name.$key is too long"))
    return String(v)
end
function label_position(o, name, default)
    p = get(o, :position, default)
    p in LABEL_POSITIONS || throw(ArgumentError("$name.position must be one of $(join(LABEL_POSITIONS, ", "))"))
    return p
end

"""
    experiment_label(spec; lumi="", annotation="")

Resolve the spec's `cms_label`/`atlas_label` (and the CMS default of
`style=:CMS`) into an `ExperimentLabel`, or `nothing`. `lumi` and
`annotation` are a comparison panel's free texts: without an explicit
luminosity/energy the label shows `lumi`, and `annotation` becomes a line below
the label. Labels of 2-D/3-D views go above the frame (`:outside`) unless a
position is given, so they do not cover coloured cells.
"""
function experiment_label(spec; lumi="", annotation="")
    default = spec isa Union{Graph2DSpec,Graph3DSpec} ? :outside : :inside_left
    cms, atlas = spec.cms_label, spec.atlas_label
    for (value, name) in ((cms, :cms_label), (atlas, :atlas_label))
        value === nothing || value isa Bool || value isa NamedTuple ||
            throw(ArgumentError("$name must be nothing, true, false, or a NamedTuple"))
    end
    explicit_cms = cms === true || cms isa NamedTuple
    explicit_atlas = atlas === true || atlas isa NamedTuple
    explicit_cms && explicit_atlas && throw(ArgumentError("choose either cms_label or atlas_label, not both"))
    # style=:CMS draws the CMS label unless it is switched off or an ATLAS label is requested.
    if explicit_cms || (cms === nothing && !explicit_atlas && style_is(spec.style, "CMS"))
        o = label_options(cms === nothing ? true : cms, CMS_LABEL_KEYS, :cms_label)
        extra = label_string(o, :extra, "Preliminary", :cms_label)
        l = label_string(o, :lumi, "", :cms_label); e = label_string(o, :energy, "", :cms_label)
        right = !isempty(l) && !isempty(e) ? "$l ($e)" : !isempty(l) ? l : !isempty(e) ? "($e)" : String(lumi)
        lines = isempty(annotation) ? String[] : [String(annotation)]
        return ExperimentLabel("CMS", label_position(o, :cms_label, default), extra, right, lines)
    elseif explicit_atlas
        o = label_options(atlas, ATLAS_LABEL_KEYS, :atlas_label)
        extra = label_string(o, :text, "Internal", :atlas_label)
        e = label_string(o, :energy, "", :atlas_label); l = label_string(o, :lumi, "", :atlas_label)
        line = join(filter(!isempty, [isempty(e) ? "" : "#sqrt{s} = $e", l]), ", ")
        isempty(line) && (line = String(lumi))
        lines = filter(!isempty, [line, String(annotation)])
        return ExperimentLabel("ATLAS", label_position(o, :atlas_label, default), extra, "", lines)
    end
    return nothing
end

"""
Label fields of a 2-D/3-D spec drawn through a converted GraphSpec scene (figure
panels): the 2-D/3-D default position (`:outside`) is made explicit.
"""
function nd_label_fields(spec)
    outside(v) = v === true ? (position=:outside,) :
                 v isa NamedTuple && !haskey(v, :position) ? merge(v, (position=:outside,)) : v
    cms = spec.cms_label === nothing && spec.atlas_label in (nothing, false) && style_is(spec.style, "CMS") ?
          (position=:outside,) : outside(spec.cms_label)
    return (cms_label=cms, atlas_label=outside(spec.atlas_label))
end

has_style_or_label(spec) = spec.style !== nothing || experiment_label(spec) !== nothing

"""Validate a spec's `style`, `cms_label` and `atlas_label` fields."""
function validate_style_and_label(spec)
    check_style_name(spec.style)
    experiment_label(spec)
    return spec
end

label_field(text) = isempty(text) ? "-" : scene_text(text)
function write_experiment_label(io, label)
    label === nothing && return
    println(io, join(("EXPLABEL", label.kind, label.position, label_field(label.extra), label_field(label.right),
                      length(label.lines), (scene_text(x) for x in label.lines)...), '\t'))
end

"""Copy of a keyword-constructed spec with some fields replaced."""
with_spec_fields(spec::T; kwargs...) where {T} =
    T(; merge(NamedTuple{fieldnames(T)}(Tuple(getfield(spec, f) for f in fieldnames(T))), (; kwargs...))...)

"""
    root_styles(; root_config=nothing)

Names accepted by the `style` option: the `TStyle`s in ROOT's style list
(`gROOT->GetListOfStyles()`, e.g. "ATLAS", "Modern", "Plain") plus the
package's "CMS" style, as reported by the installed native ROOT.
"""
function root_styles(; root_config=nothing)
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    return String.(split(strip(read(Cmd([exe, "--styles"]), String)), '\n'))
end
