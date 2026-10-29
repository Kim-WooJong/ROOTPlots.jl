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

"""Caller-provided panel labels. `annotation` and `lumi` are descriptive text, never inferred physics claims; object owns only text."""
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
end

function PanelSpec(xlabel::AbstractString; ylabel="Events", annotation="", lumi="",
                   logx=false, logy=false, ratio_logy=false, ymin=nothing, ymax=nothing,
                   ratio_limits=nothing)
    limits = ratio_limits === nothing ? nothing : (Float64(ratio_limits[1]), Float64(ratio_limits[2]))
    return PanelSpec(String(xlabel), String(ylabel), String(annotation), String(lumi),
                     logx, logy, ratio_logy, isnothing(ymin) ? nothing : Float64(ymin),
                     isnothing(ymax) ? nothing : Float64(ymax), limits)
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
