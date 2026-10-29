"""Resolved control-panel limits. Log limits are physical coordinates, never logarithms."""
function panel_limits(cmp::Comparison, spec::PanelSpec, fallback_ratio)
    edges = cmp.data.edges
    spec.logx && first(edges) <= 0 &&
        throw(ArgumentError("logx requires every displayed bin edge to be positive"))
    candidate = spec.ratio_limits === nothing ? fallback_ratio : spec.ratio_limits
    length(candidate) == 2 || throw(ArgumentError("ratio_limits must contain two bounds"))
    rlo, rhi = Float64.(candidate)
    isfinite(rlo) && isfinite(rhi) && rlo < rhi && (!spec.ratio_logy || rlo > 0) ||
        throw(ArgumentError("invalid ratio limits for selected scale"))
    upper = max(maximum(cmp.prediction .+ cmp.prediction_errors.high),
                maximum(cmp.data.sumw .+ cmp.data_errors.high))
    if spec.logy
        positives = Float64[]
        for values in cmp.cumulative
            append!(positives, filter(x -> isfinite(x) && x > 0, values))
        end
        for i in eachindex(cmp.prediction)
            for x in (cmp.prediction[i], cmp.prediction[i] - cmp.prediction_errors.low[i],
                      cmp.prediction[i] + cmp.prediction_errors.high[i],
                      cmp.data.sumw[i], cmp.data.sumw[i] - cmp.data_errors.low[i],
                      cmp.data.sumw[i] + cmp.data_errors.high[i])
                isfinite(x) && x > 0 && push!(positives, x)
            end
        end
        isempty(positives) && (spec.ymin === nothing || spec.ymax === nothing) &&
            throw(ArgumentError("automatic logy limits require positive drawable values"))
        lo = spec.ymin === nothing ? minimum(positives) / 2 : spec.ymin
        hi = spec.ymax === nothing ? max(upper, maximum(positives)) * 1.46 : spec.ymax
    else
        lo = spec.ymin === nothing ? 0.0 : spec.ymin
        hi = spec.ymax === nothing ? max(1.0, upper * 1.46) : spec.ymax
    end
    isfinite(lo) && isfinite(hi) && lo < hi && (!spec.logy || lo > 0) ||
        throw(ArgumentError("invalid y limits for selected scale"))
    return (; ymin=lo, ymax=hi, ratio_min=rlo, ratio_max=rhi)
end

function validate_style(style::PlotStyle)
    for color in (style.marker_color, style.band_color)
        color === nothing || occursin(r"^#[0-9a-fA-F]{6}$", color) ||
            throw(ArgumentError("style colors must use #RRGGBB"))
    end
    if style.marker_style !== nothing
        1 <= style.marker_style <= 50 || throw(ArgumentError("invalid ROOT marker style"))
    end
    style.marker_size === nothing ||
        (isfinite(style.marker_size) && style.marker_size > 0) ||
        throw(ArgumentError("marker_size must be positive and finite"))
    style.line_width === nothing || style.line_width > 0 ||
        throw(ArgumentError("line_width must be positive"))
    for code in (style.stack_fill_style, style.band_fill_style)
        code === nothing || code in (0, 1001, 3345) ||
            throw(ArgumentError("supported fill styles are 0, 1001, and 3345"))
    end
    style.legend_columns === nothing || style.legend_columns > 0 ||
        throw(ArgumentError("legend_columns must be positive"))
    return style
end

"""Validate output selection and return stable PNG/PDF keys, with `nothing` for unselected formats."""
function plot_output_paths(outstem::AbstractString, format)
    format in (:png, :pdf, :both) ||
        throw(ArgumentError("format must be :png, :pdf, or :both"))
    return (; png=format in (:png, :both) ? outstem * ".png" : nothing,
              pdf=format in (:pdf, :both) ? outstem * ".pdf" : nothing)
end

"""Move only requested native output files; the temporary directory discards the other format."""
function move_plot_outputs(stem::AbstractString, paths)
    for (extension, destination) in ((".png", paths.png), (".pdf", paths.pdf))
        destination === nothing && continue
        source = stem * extension
        isfile(source) || error("native ROOT renderer did not create $extension")
        mv(source, destination; force=true)
    end
    return paths
end
