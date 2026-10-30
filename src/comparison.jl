"""Absolute lower/upper bin error magnitudes, not interval endpoints. Inputs are nonnegative finite vectors; no external state is owned."""
struct BinErrors
    low::Vector{Float64}
    high::Vector{Float64}

    function BinErrors(low::AbstractVector{<:Real}, high::AbstractVector{<:Real})
        lo, hi = Float64.(low), Float64.(high)
        length(lo) == length(hi) && all(isfinite, lo) && all(isfinite, hi) &&
            all(lo .>= 0) && all(hi .>= 0) || throw(ArgumentError("invalid bin errors"))
        new(lo, hi)
    end
end

BinErrors(sym::AbstractVector{<:Real}) = BinErrors(sym, sym)

"""One nonnegative stack component. Input histogram's signed event weights have already been summed; visible negative yields cannot be stacked conventionally."""
struct Component
    label::String
    hist::Hist1D
    color::String

    function Component(label::AbstractString, hist::Hist1D, color::AbstractString)
        isempty(label) && throw(ArgumentError("empty component label"))
        occursin(r"^#[0-9a-fA-F]{6}$", color) || throw(ArgumentError("use a #RRGGBB color"))
        all(hist.sumw .>= 0) || throw(ArgumentError("negative component-bin yield cannot be stacked"))
        new(String(label), hist, String(color))
    end
end

"""
    PanelSpec(xlabel; ylabel="Events", annotation="", lumi="", logx=false, logy=false,
              ratio_logy=false, ymin=nothing, ymax=nothing, ratio_limits=nothing,
              xaxis=AxisStyle(), yaxis=AxisStyle(), ratio_yaxis=AxisStyle(), canvas=CanvasStyle(),
              style=nothing, cms_label=nothing, atlas_label=nothing)

One stacked comparison panel with its ratio pad. Caller-provided labels:
`annotation` and `lumi` are descriptive text, never inferred physics claims.

- `xlabel`, `ylabel`: axis titles (x on the ratio pad, y on the upper pad).
- `logx` (both pads), `logy` (upper pad), `ratio_logy` (ratio pad): `TPad::SetLogx/SetLogy`.
- `ymin`, `ymax`, `ratio_limits`: physical-coordinate ranges (never logarithms).
- `xaxis`: [`AxisStyle`](@ref) of the shared x axis (drawn on the ratio pad; the upper
  pad keeps hidden labels); `yaxis`: upper-pad y axis; `ratio_yaxis`: ratio-pad y axis.
  `xaxis.bin_labels` label the comparison bins.
- `canvas`: [`CanvasStyle`](@ref) with `grid_x`, `grid_y`, `ticks_x`, `ticks_y` (both
  pads). The comparison layout owns margins, titles, and global settings.
- `style`, `cms_label`, `atlas_label`: ROOT style and experiment label as in
  [`GraphSpec`](@ref) (panels of one figure must name the same style). With a style
  or a label, both pads take the style's margins, ticks and text sizes rescaled
  to the pad heights, and the legend moves into the frame. In every case the
  ratio pad's text has the upper y axis's pixel size and the two y titles line
  up; explicit `AxisStyle` sizes are ROOT's values for their own pad (see the API
  documentation, "Two-pad views"). A label carries the
  panel texts: `lumi` is the CMS luminosity text (or the ATLAS second line) unless
  the label gives `lumi`/`energy`, and `annotation` is drawn as a line below it.
"""
struct PanelSpec
    xlabel::String
    ylabel::String
    annotation::String
    lumi::String
    logx::Bool
    logy::Bool
    ratio_logy::Bool
    ymin::Union{Nothing,Float64}
    ymax::Union{Nothing,Float64}
    ratio_limits::Union{Nothing,Tuple{Float64,Float64}}
    xaxis::AxisStyle
    yaxis::AxisStyle
    ratio_yaxis::AxisStyle
    canvas::CanvasStyle
    style::Union{Nothing,Symbol}
    cms_label::Union{Nothing,Bool,NamedTuple}
    atlas_label::Union{Nothing,Bool,NamedTuple}
end

# The original ten- and fourteen-field positional forms remain valid.
PanelSpec(xlabel, ylabel, annotation, lumi, logx, logy, ratio_logy, ymin, ymax, ratio_limits) =
    PanelSpec(xlabel, ylabel, annotation, lumi, logx, logy, ratio_logy, ymin, ymax, ratio_limits,
              AxisStyle(), AxisStyle(), AxisStyle(), CanvasStyle())
PanelSpec(xlabel, ylabel, annotation, lumi, logx, logy, ratio_logy, ymin, ymax, ratio_limits,
          xaxis, yaxis, ratio_yaxis, canvas) =
    PanelSpec(xlabel, ylabel, annotation, lumi, logx, logy, ratio_logy, ymin, ymax, ratio_limits,
              xaxis, yaxis, ratio_yaxis, canvas, nothing, nothing, nothing)

function PanelSpec(xlabel::AbstractString; ylabel="Events", annotation="", lumi="",
                   logx=false, logy=false, ratio_logy=false, ymin=nothing, ymax=nothing,
                   ratio_limits=nothing, xaxis=AxisStyle(), yaxis=AxisStyle(),
                   ratio_yaxis=AxisStyle(), canvas=CanvasStyle(), style=nothing, cms_label=nothing,
                   atlas_label=nothing)
    limits = ratio_limits === nothing ? nothing : (Float64(ratio_limits[1]), Float64(ratio_limits[2]))
    return PanelSpec(String(xlabel), String(ylabel), String(annotation), String(lumi),
                     logx, logy, ratio_logy, isnothing(ymin) ? nothing : Float64(ymin),
                     isnothing(ymax) ? nothing : Float64(ymax), limits, xaxis, yaxis, ratio_yaxis, canvas,
                     style, cms_label, atlas_label)
end

const PANEL_CANVAS_FIELDS = (:grid_x, :grid_y, :ticks_x, :ticks_y)
has_panel_axes(spec::PanelSpec) = spec.xaxis != AxisStyle() || spec.yaxis != AxisStyle() ||
    spec.ratio_yaxis != AxisStyle() || spec.canvas != CanvasStyle()

"""Validate the comparison-panel axis styles for a panel with `nbins` bins."""
function validate_panel_axes(spec::PanelSpec, nbins)
    validate_axis(spec.xaxis; name=:x); validate_axis(spec.yaxis; name=:y)
    validate_axis(spec.ratio_yaxis; name=:ratio)
    validate_canvas(spec.canvas)
    unsupported = [f for f in fieldnames(CanvasStyle) if !(f in PANEL_CANVAS_FIELDS) && getfield(spec.canvas, f) !== nothing]
    isempty(unsupported) ||
        throw(ArgumentError("PanelSpec canvas supports $(join(PANEL_CANVAS_FIELDS, ", ")); the comparison layout owns $(join(unsupported, ", "))"))
    validate_style_and_label(spec)
    check_axis_bins(spec.xaxis, nbins, :x, "the comparison panel")
    check_axis_bins(spec.yaxis, nothing, :y, "the comparison panel")
    check_axis_bins(spec.ratio_yaxis, nothing, :ratio, "the comparison panel")
    ys = filter(!isnothing, [spec.yaxis.exponent_offset, spec.ratio_yaxis.exponent_offset])
    length(unique(ys)) <= 1 ||
        throw(ArgumentError("yaxis and ratio_yaxis exponent_offset differ (TGaxis::SetExponentOffset is global)"))
    return spec
end

PanelSpec(xlabel::AbstractString, ylabel::AbstractString, annotation::AbstractString,
          lumi::AbstractString) = PanelSpec(xlabel; ylabel, annotation, lumi)

"""Native ROOT marker/fill codes; `nothing` preserves the renderer's defaults."""
Base.@kwdef struct PlotStyle
    marker_style::Union{Nothing,Int} = nothing
    marker_size::Union{Nothing,Float64} = nothing
    marker_color::Union{Nothing,String} = nothing
    line_width::Union{Nothing,Int} = nothing
    stack_fill_style::Union{Nothing,Int} = nothing
    band_fill_style::Union{Nothing,Int} = nothing
    band_color::Union{Nothing,String} = nothing
    legend_columns::Union{Nothing,Int} = nothing
end

"""Numerical stack/data comparison; ratio slots are `NaN` only where prediction ≤0, with `valid` identifying drawable bins. No rendering or I/O."""
struct Comparison
    data::Hist1D
    components::Vector{Component}
    data_errors::BinErrors
    prediction_errors::BinErrors
    prediction::Vector{Float64}
    cumulative::Vector{Vector{Float64}}
    ratio::Vector{Float64}
    ratio_data_errors::BinErrors
    ratio_prediction_errors::BinErrors
    valid::BitVector
    uncertainty_label::String
end

"""Compute stack totals, independent statistical default uncertainty, and Data/prediction ratio. Supplied prediction errors replace that default; no systematics are invented or double-counted."""
function comparison(data::Hist1D, components::AbstractVector{Component};
                    data_errors::BinErrors=BinErrors(sqrt.(data.sumw2)),
                    prediction_errors::Union{Nothing,BinErrors}=nothing,
                    uncertainty_label::AbstractString=prediction_errors === nothing ? "Prediction stat. unc." : "Prediction uncertainty")
    isempty(components) && throw(ArgumentError("at least one prediction component is needed"))
    all(c.hist.edges == data.edges for c in components) || throw(ArgumentError("bin edges disagree"))
    all(data.sumw .>= 0) || throw(ArgumentError("observed data bin yields must be nonnegative"))
    n = length(data.sumw)
    length(data_errors.low) == n || throw(ArgumentError("data error length disagrees"))
    pred = zeros(n)
    statvar = zeros(n)
    cumulative = Vector{Vector{Float64}}()
    for c in components
        pred .+= c.hist.sumw
        statvar .+= c.hist.sumw2
        push!(cumulative, copy(pred))
    end
    pe = prediction_errors === nothing ? BinErrors(sqrt.(statvar)) : prediction_errors
    length(pe.low) == n || throw(ArgumentError("prediction error length disagrees"))
    valid = BitVector(pred .> 0)
    ratio = fill(NaN, n)
    dlo, dhi, plo, phi = (zeros(n) for _ in 1:4)
    for i in eachindex(pred)
        valid[i] || continue
        ratio[i] = data.sumw[i] / pred[i]
        dlo[i], dhi[i] = data_errors.low[i] / pred[i], data_errors.high[i] / pred[i]
        plo[i], phi[i] = pe.low[i] / pred[i], pe.high[i] / pred[i]
    end
    return Comparison(data, collect(components), data_errors, pe, pred, cumulative,
                      ratio, BinErrors(dlo, dhi), BinErrors(plo, phi),
                      valid, String(uncertainty_label))
end
