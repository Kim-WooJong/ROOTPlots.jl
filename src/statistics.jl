"""Validate the visible-bin arrays again because `Hist1D` owns mutable vectors."""
function validate_stats_hist(h::Hist1D)
    n = length(h.sumw)
    n > 0 && length(h.edges) == n + 1 && length(h.sumw2) == n &&
        all(isfinite, h.edges) && all(diff(h.edges) .> 0) &&
        all(isfinite, h.sumw) && all(isfinite, h.sumw2) && all(>=(0), h.sumw2) ||
        throw(ArgumentError("invalid histogram arrays"))
    return n
end

"""Merge complete visible bins onto an exact subset of the original edges.

The outer edges must be unchanged: no fractional-bin splitting, cropping, or
implicit movement of visible content into flow bins is performed. Independent
per-bin `sumw2` values add in quadrature; signed `sumw` values are allowed.
"""
function rebin(h::Hist1D, new_edges::AbstractVector{<:Real})
    validate_stats_hist(h)
    edges = Float64.(new_edges)
    length(edges) >= 2 && all(isfinite, edges) && all(diff(edges) .> 0) ||
        throw(ArgumentError("new bin edges must be finite and strictly increasing"))
    first(edges) == first(h.edges) && last(edges) == last(h.edges) ||
        throw(ArgumentError("rebinning must preserve the outer histogram edges"))
    positions = searchsortedfirst.(Ref(h.edges), edges)
    all(i <= length(h.edges) && h.edges[i] == edge for (i, edge) in zip(positions, edges)) ||
        throw(ArgumentError("new edges must exactly match existing bin edges"))
    values = [sum(@view h.sumw[positions[j]:positions[j+1]-1])
              for j in 1:length(edges)-1]
    variances = [sum(@view h.sumw2[positions[j]:positions[j+1]-1])
                 for j in 1:length(edges)-1]
    return Hist1D(edges, values, variances)
end

"""Visible-bin sum and propagated error. `width=true` integrates bin contents as densities."""
function integral(h::Hist1D; bins=eachindex(h.sumw), width::Bool=false)
    n = validate_stats_hist(h)
    indices = collect(bins)
    all(i -> i isa Integer && 1 <= i <= n, indices) &&
        all(diff(indices) .> 0) || throw(ArgumentError("bins must be increasing valid indices"))
    value = sum(i -> h.sumw[i] * (width ? h.edges[i+1]-h.edges[i] : 1.0), indices;
                init=0.0)
    variance = sum(i -> h.sumw2[i] * (width ? (h.edges[i+1]-h.edges[i])^2 : 1.0),
                   indices; init=0.0)
    isfinite(value) && isfinite(variance) || throw(ArgumentError("integral overflows Float64"))
    return (value=value, error=sqrt(variance))
end

"""ROOT-style effective entries `(Σw)^2/Σw²` for visible bins, not the fill count."""
function effective_entries(h::Hist1D)
    validate_stats_hist(h)
    total = sum(h.sumw)
    variance = sum(h.sumw2)
    isfinite(total) && isfinite(variance) || throw(ArgumentError("histogram total overflows Float64"))
    result = variance == 0 ? abs(total) : (abs(total) / sqrt(variance))^2
    isfinite(result) || throw(ArgumentError("effective entries overflow Float64"))
    return result
end

function binned_average(h::Hist1D)
    n = validate_stats_hist(h)
    total = sum(h.sumw)
    isfinite(total) && total > 0 || throw(ArgumentError("binned moments need positive Σw"))
    centers = [(h.edges[i] / 2 + h.edges[i+1] / 2) for i in 1:n]
    all(isfinite, centers) || throw(ArgumentError("bin centers overflow Float64"))
    mean = sum(h.sumw .* centers) / total
    isfinite(mean) || throw(ArgumentError("binned mean overflows Float64"))
    return (mean=mean, total=total, centers=centers)
end

function binned_moments(h::Hist1D)
    average = binned_average(h)
    (; mean, total, centers) = average
    variance = sum(h.sumw .* ((centers .- mean) .^ 2)) / total
    isfinite(mean) && isfinite(variance) && variance >= 0 ||
        throw(ArgumentError("binned variance is negative or nonfinite (signed weights?)"))
    return (mean=mean, std=sqrt(variance))
end

"""Mean of visible bin centres weighted by `Σw`; not the original-event mean."""
binned_mean(h::Hist1D) = binned_average(h).mean

"""Population standard deviation of visible bin centres; no Bessel or bin-width correction."""
binned_std(h::Hist1D) = binned_moments(h).std

"""ROOT-style approximate mean error, `binned_std(h)/√effective_entries(h)`."""
function binned_mean_error(h::Hist1D)
    n = effective_entries(h)
    n > 0 || throw(ArgumentError("mean error needs positive effective entries"))
    return binned_std(h) / sqrt(n)
end

"""One visible-bin summary. Flow `(Σw, Σw²)` pairs are reported separately, not included."""
function histstats(h::Hist1D; underflow=nothing, overflow=nothing)
    for flow in (underflow, overflow)
        flow === nothing && continue
        flow isa Tuple && length(flow) == 2 &&
            all(x -> x isa Real && isfinite(x), flow) && flow[2] >= 0 ||
            throw(ArgumentError("flow must be a finite (Σw, Σw²) pair"))
    end
    moments = binned_moments(h)
    neff = effective_entries(h)
    return (integral=integral(h), effective_entries=neff,
            mean=moments.mean, mean_error=moments.std/sqrt(neff), std=moments.std,
            negative_bins=any(<(0), h.sumw), underflow=underflow, overflow=overflow)
end

"""Parameter correlation matrix. Fixed/zero-error parameter rows and columns are `NaN`."""
function fit_correlation(fit::FitResult)
    n = length(fit.parameters)
    size(fit.covariance) == (n,n) && length(fit.errors) == n && length(fit.fixed) == n ||
        throw(ArgumentError("invalid fit covariance dimensions"))
    result = fill(NaN, n, n)
    for i in 1:n, j in 1:n
        if !fit.fixed[i] && !fit.fixed[j] && fit.errors[i] > 0 && fit.errors[j] > 0
            result[i,j] = fit.covariance[i,j] / (fit.errors[i]*fit.errors[j])
        end
    end
    return result
end

"""Observed-minus-model bin-centre residuals and pulls for a χ² histogram fit.

Bins outside the fit range or with zero error have `valid=false` and `NaN` residual/pull.
"""
function fit_residuals(h::Hist1D, fit::FitResult)
    n = validate_stats_hist(h)
    fit.chi2_kind == :chi2 || throw(ArgumentError("pulls require a χ² histogram fit"))
    centers = [(h.edges[i] / 2 + h.edges[i+1] / 2) for i in 1:n]
    predicted = [model_value(fit.model, fit.parameters, x) for x in centers]
    all(isfinite, predicted) || throw(ArgumentError("fit prediction is nonfinite"))
    valid = [fit.range[1] <= centers[i] <= fit.range[2] && h.sumw2[i] > 0 for i in 1:n]
    residual = [valid[i] ? h.sumw[i]-predicted[i] : NaN for i in 1:n]
    pull = [valid[i] ? residual[i]/sqrt(h.sumw2[i]) : NaN for i in 1:n]
    return (centers=centers, observed=copy(h.sumw), predicted=predicted,
            residual=residual, pull=pull, valid=valid)
end

"""Bin-centre points with `sqrt(sumw2)` vertical errors for a graph overlay.

This is not a step histogram. Bin half-widths are deliberately not supplied as
x errors: they are not measurement uncertainties, and ROOT fits use x errors.
Fit the original `Hist1D`, not this display-only graph, to retain TH1 semantics.
"""
function histogram_points(h::Hist1D; label="Histogram", style=GraphStyle())
    n = validate_stats_hist(h)
    centers = [h.edges[i] / 2 + h.edges[i+1] / 2 for i in 1:n]
    all(isfinite, centers) || throw(ArgumentError("bin centers overflow Float64"))
    return GraphSeries(centers, h.sumw; yerrors=sqrt.(h.sumw2), label, style)
end

"""Display visible bin contents as a step line with exact (possibly variable) bin edges.

This is drawing geometry, not a new histogram or fit input. It carries no
statistical errors; overlay `histogram_points(h)` when error bars are wanted.
"""
function histogram_steps(h::Hist1D; label="Histogram", style=GraphStyle())
    n = validate_stats_hist(h)
    x = Vector{Float64}(undef, 2n)
    y = Vector{Float64}(undef, 2n)
    for i in 1:n
        x[2i-1] = h.edges[i]
        x[2i] = h.edges[i+1]
        y[2i-1] = y[2i] = h.sumw[i]
    end
    return GraphSeries(x, y; label, style, draw=:line, errors=:none)
end

"""Render a standalone step histogram with a visible zero baseline on linear y axes."""
function render_root(h::Hist1D, spec::GraphSpec, outstem::AbstractString;
                     draw=:steps, label="Histogram", style=GraphStyle(), format=:png,
                     stats=false, attached_fit=nothing, root_config=nothing)
    if draw isa AbstractString && uppercase(strip(draw)) != "STEPS"
        return render_hist_root(h, spec, outstem;
                                mode=String(draw), label, style, format, stats,
                                attached_fit, root_config)
    end
    stats && throw(ArgumentError("statistics boxes require native ROOT histogram drawing"))
    attached_fit === nothing ||
        throw(ArgumentError("attached_fit requires a raw ROOT draw-option string"))
    mode = hist_draw_mode(draw)
    mode == :steps || return render_hist_root(h, spec, outstem;
                                            mode, label, style, format, root_config)
    line = histogram_steps(h; label, style)
    ylimits = spec.ylimits
    if ylimits === nothing && !spec.logy
        lo, hi = min(0.0, minimum(h.sumw)), max(0.0, maximum(h.sumw))
        if lo == hi
            ylimits = (-0.5, 0.5)
        else
            span = hi - lo
            isfinite(span) || throw(ArgumentError("histogram y range overflows Float64"))
            pad = 0.05span
            ylimits = (lo - pad, hi + pad)
        end
    end
    display_spec = with_graph_spec(spec;
                             xlimits=spec.xlimits === nothing ?
                                     (first(h.edges), last(h.edges)) : spec.xlimits,
                             ylimits)
    return render_root([line], display_spec, outstem; format, root_config)
end

"""Display valid χ² histogram-fit pulls as graph points; invalid bins are omitted."""
function fit_pull_points(h::Hist1D, fit::FitResult; label="Pull", style=GraphStyle())
    fit.status == 0 && fit.valid || throw(ArgumentError("pull plot needs a valid fit"))
    residuals = fit_residuals(h, fit)
    indices = findall(residuals.valid)
    isempty(indices) && throw(ArgumentError("fit has no valid pull bins"))
    return GraphSeries(residuals.centers[indices], residuals.pull[indices];
                       label, style, errors=:none)
end

"""Render a χ² histogram fit and its valid-bin pulls in aligned ROOT pads.

The upper pad shows every visible bin, while the curve and pull points use the
fit's range. The caller must supply the result of `fit_root(h, ...)`; FitResult
does not retain source-data identity. `spec.logy` applies only to the upper pad.
"""
function render_root(h::Hist1D, fit::FitResult, spec::GraphSpec, outstem::AbstractString;
                     pull_limits=nothing, curve_samples::Int=300,
                     data_label="Histogram", fit_label="ROOT fit",
                     data_style=GraphStyle(), fit_style=GraphStyle(color="#cc3311"),
                     pull_style=GraphStyle(), format=:png, root_config=nothing)
    paths = plot_output_paths(outstem, format)
    has_scene(spec) && throw(ArgumentError("scene styling is currently supported on one-pad plots, not fit/pull layouts"))
    validate_stats_hist(h)
    all(x -> 0 < x <= 10000, spec.size) ||
        throw(ArgumentError("canvas dimensions must be in 1:10000"))
    pull = fit_pull_points(h, fit; style=pull_style)
    curve = fit_curve(fit; n=curve_samples, label=fit_label, style=fit_style)
    data = histogram_points(h; label=data_label, style=data_style)
    xlimits = spec.xlimits === nothing ? (first(h.edges), last(h.edges)) : spec.xlimits
    display_spec = GraphSpec(; xlabel=spec.xlabel, ylabel=spec.ylabel, title=spec.title,
                             xlimits, ylimits=spec.ylimits, logx=spec.logx,
                             logy=spec.logy, legend=spec.legend, size=spec.size)
    limits = graph_limits([curve, data], display_spec)
    span = max(3.0, 1.2 * maximum(abs, pull.y))
    ylim = pull_limits === nothing ? (-span, span) : Tuple(Float64.(pull_limits))
    length(ylim) == 2 && all(isfinite, ylim) && ylim[1] < 0 < ylim[2] ||
        throw(ArgumentError("pull limits must be finite, increasing, and contain zero"))
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        input = joinpath(dir,"fit-plot.tsv")
        write_fit_plot_spec(input,[curve,data,pull],display_spec,limits,ylim)
        stem = joinpath(dir,"rendered")
        run(Cmd([exe,input,stem]))
        move_plot_outputs(stem,paths)
    end
    return paths
end

"""Asymptotic upper-tail χ² probability for valid χ² or unweighted Poisson fits.

Assumes suitable bin counts, a preselected fit range, no active parameter bound,
and a model sufficiently regular for the χ² approximation. Weighted Poisson and
fits with ROOT-normalized errors are deliberately refused.
"""
function chi2_pvalue(fit::FitResult)
    fit.status == 0 && fit.valid && fit.ndf > 0 &&
        isfinite(fit.chi2) && fit.chi2 >= 0 && !fit.normalized_errors &&
        fit.chi2_kind in (:chi2, :poisson_deviance) ||
        throw(ArgumentError("χ² p-value is not defined for this fit result"))
    return gamma_inc(fit.ndf/2, fit.chi2/2)[2]
end
