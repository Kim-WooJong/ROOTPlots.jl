# Safe numeric operations on Hist1D/Hist2D/Hist3D: scalar scaling, merging of
# independent samples, and projections. Julia owns these numbers; ROOT only draws
# the results. Every operation returns new arrays and never mutates its input.

const AnyHist = Union{Hist1D,Hist2D,Hist3D}

# ---------------------------------------------------------------------------
# Full arrays (visible + flow bins) and reconstruction

"""Revalidate a histogram's mutable arrays, including the stored flow arrays."""
function validate_hist_arrays(h::Hist1D)
    validate_stats_hist(h)
    return h
end

function validate_hist_arrays(h::Union{Hist2D,Hist3D})
    validate_stats_hist_nd(h)
    f = h.bins_with_flow
    if f !== nothing
        edges = hist_nd_edges(h)
        shape = map(e -> length(e) + 1, edges)
        visible = map(e -> 2:length(e), edges)
        size(f.sumw) == shape && size(f.sumw2) == shape && all(isfinite, f.sumw) &&
            all(isfinite, f.sumw2) && all(>=(0), f.sumw2) &&
            f.sumw[visible...] == h.sumw && f.sumw2[visible...] == h.sumw2 ||
            throw(ArgumentError("histogram arrays were modified: bins_with_flow no longer matches sumw/sumw2"))
    end
    return h
end

"""New `(sumw, sumw2)` arrays with flow bins (index 1 and `end` on every axis).

Flow bins are zero when the histogram stores none (`bins_with_flow === nothing`).
"""
function full_arrays(h::Hist1D)
    return ([h.underflow[1]; h.sumw; h.overflow[1]], [h.underflow[2]; h.sumw2; h.overflow[2]])
end

function full_arrays(h::Union{Hist2D,Hist3D})
    h.bins_with_flow === nothing || return (copy(h.bins_with_flow.sumw), copy(h.bins_with_flow.sumw2))
    edges = hist_nd_edges(h)
    shape = map(e -> length(e) + 1, edges)
    visible = map(e -> 2:length(e), edges)
    w, v = zeros(shape), zeros(shape)
    w[visible...] .= h.sumw
    v[visible...] .= h.sumw2
    return w, v
end

hist_edges(h::Hist1D) = (h.edges,)
hist_edges(h::Union{Hist2D,Hist3D}) = hist_nd_edges(h)

"""Build a histogram of `h`'s type from full arrays; the visible part is taken from them."""
function from_full(h::Hist1D, w, v, keep_flow::Bool, entries, moments)
    return Hist1D(h.edges, w[2:end-1], v[2:end-1]; underflow=(w[1], v[1]),
                  overflow=(w[end], v[end]), entries, moments)
end

function from_full(h::Union{Hist2D,Hist3D}, w, v, keep_flow::Bool, entries, moments)
    visible = map(e -> 2:length(e), hist_nd_edges(h))
    flow = keep_flow ? (sumw=w, sumw2=v) : nothing
    return typeof(h)(hist_nd_edges(h)..., w[visible...], v[visible...];
                     bins_with_flow=flow, entries, moments)
end

has_flow_arrays(h::Hist1D) = true
has_flow_arrays(h::Union{Hist2D,Hist3D}) = h.bins_with_flow !== nothing

# ---------------------------------------------------------------------------
# Scalar scaling (TH1::Scale)

function scale_factor(c::Real)
    f = Float64(c)
    isfinite(f) || throw(ArgumentError("histogram scale factor must be finite, got $c"))
    f2 = f * f
    isfinite(f2) || throw(ArgumentError("histogram scale factor $c squared overflows Float64"))
    return f, f2
end

"""Multiply every element by `c`; normalize -0.0 to 0.0; error on Float64 overflow."""
function scaled_values(a::AbstractArray{Float64}, c::Float64, what::AbstractString)
    out = similar(a)
    for i in eachindex(a)
        y = a[i] * c
        isfinite(y) || throw(ArgumentError("scaling $what overflows Float64"))
        out[i] = y == 0 ? 0.0 : y
    end
    return out
end

function scaled_moments(m, c::Float64, c2::Float64)
    m === nothing && return nothing
    # TH1::Scale: index 2 (Σw²) scales with c², every other statistic with c.
    s = scaled_values(collect(m), c, "the histogram moments")
    s[2] = only(scaled_values([m[2]], c2, "the histogram moments"))
    return Tuple(s)
end

"""
    c * h, h * c, h / c

Scale a `Hist1D`, `Hist2D`, or `Hist3D` by a finite real number, like ROOT's
`TH1::Scale(c)`: every bin including under/overflow gets `Σw → c Σw` and
`Σw² → c² Σw²`; the unbinned moments (`TH1::GetStats` order) are scaled the same
way (`Σw²` by `c²`, all others by `c`), so the mean and standard deviation are
unchanged; `entries` (the fill count) is unchanged. A new histogram with its own
arrays is returned; the input is not modified.

This equals a direct fill of the same events with every weight multiplied by
`c` (up to rounding), and is the operation for reweighting one sample (cross
section × luminosity / generated events, a luminosity scale, ...).

- `c < 0` is allowed (ROOT keeps the statistics, as here).
- `c == 0` gives all-zero contents, variances, and moments with the fill count
  kept, as `Scale(0)` does in ROOT (ROOT then shows Mean = Std Dev = 0).
- `h / c` is `h * inv(c)` (ROOT `Scale(1/c)`); `c == 0` throws `ArgumentError`.
- Non-finite factors and results that overflow Float64 throw `ArgumentError`.
- Ordinary Float64 limits apply at the small end: `c²` and the products are
  not checked for underflow, so a tiny factor (`c²` is subnormal below
  `abs(c) ≈ 1.5e-154` and zero below `≈ 1.5e-162`) or tiny contents can round
  variances and moments to reduced-precision values or zero while `Σw` stays
  nonzero, as in ROOT's own `Scale`.

There is deliberately no normalization or density helper: scaling by a
data-dependent number such as `1/integral(h)` correlates the bins, which
per-bin `Σw²` cannot describe.
"""
function Base.:*(c::Real, h::AnyHist)
    validate_hist_arrays(h)
    f, f2 = scale_factor(c)
    w, v = full_arrays(h)
    sw = scaled_values(w, f, "the histogram contents")
    sv = scaled_values(v, f2, "the histogram variances")
    return from_full(h, sw, sv, has_flow_arrays(h), h.entries, scaled_moments(h.moments, f, f2))
end

Base.:*(h::AnyHist, c::Real) = c * h

function Base.:/(h::AnyHist, c::Real)
    iszero(c) && throw(ArgumentError("cannot divide a histogram by zero"))
    f = Float64(c)
    isfinite(f) || throw(ArgumentError("histogram divisor must be finite, got $c"))
    return inv(f) * h
end

# ---------------------------------------------------------------------------
# Merging independent samples (TH1::Add with c = 1 / TH1::Merge)

data_arrays(h::Hist1D) = (h.sumw, h.sumw2)
data_arrays(h::Union{Hist2D,Hist3D}) =
    h.bins_with_flow === nothing ? (h.sumw, h.sumw2) :
    (h.sumw, h.sumw2, h.bins_with_flow.sumw, h.bins_with_flow.sumw2)

shares_data(a, b) = a === b || any(x -> any(y -> x === y, data_arrays(b)), data_arrays(a))

function merged_entries(hs)
    any(h -> h.entries === nothing, hs) && return nothing
    total = 0
    for h in hs
        total, overflow = Base.Checked.add_with_overflow(total, h.entries)
        overflow && throw(ArgumentError("merged entry count overflows Int"))
    end
    return total
end

function merged_moments(hs)
    any(h -> h.moments === nothing, hs) && return nothing
    total = collect(first(hs).moments)
    for h in hs[2:end]
        total .+= h.moments
    end
    all(isfinite, total) || throw(ArgumentError("merged histogram moments overflow Float64"))
    return Tuple(total)
end

"""
    merge(a::Hist1D, b::Hist1D, more::Hist1D...)
    merge(a::Hist2D, b::Hist2D, more...), merge(a::Hist3D, b::Hist3D, more...)

Combine histograms of **independent event samples** filled with the same
binning, like ROOT's `TH1::Add(h, 1)` / `TH1::Merge`. All bins including
under/overflow add `Σw` and `Σw²`; known `entries` add; known unbinned moments
add. If any input has unknown `entries` (or unknown moments), the result's
`entries` (or moments) are unknown too, so binned and unbinned statistics are
never mixed. For 2D/3D the result keeps `bins_with_flow` when any input has
it (missing flow bins count as zero), and its visible bins are taken from that
full array. The result equals a direct fill of the concatenated events (up to
floating-point summation order). Inputs are not modified.

Checked requirements:

- same histogram type and exactly equal edges on every axis (use [`rebin`](@ref)
  first to merge onto common coarser edges);
- distinct inputs: passing the same histogram twice (or histograms sharing
  arrays) throws `ArgumentError`. To weight one sample twice, scale it
  (`2 * h`): adding a sample to itself treats the same events as independent
  and gives `2Σw²` instead of `4Σw²`.

Independence is the caller's responsibility: overlapping samples that were
copied or filled twice cannot be detected and give underestimated errors.
Signed weights are added as they are; bins can cancel to small or negative
`Σw` while `Σw²` keeps growing, and ROOT's effective entries and mean then
describe the signed sum, not an event count. There is no subtraction: ROOT's
`Add(h, -1)` resets the statistics, and a difference of correlated samples does
not have variance `Σw²₁ + Σw²₂`.
"""
function Base.merge(a::H, b::H, more::H...) where {H<:AnyHist}
    hs = (a, b, more...)
    foreach(validate_hist_arrays, hs)
    for i in eachindex(hs), j in i+1:length(hs)
        shares_data(hs[i], hs[j]) && throw(ArgumentError(
            "merge inputs $i and $j are the same histogram (or share arrays); merge adds " *
            "independent samples - to reweight one sample use scalar scaling such as 2 * h"))
    end
    edges = hist_edges(a)
    for (k, h) in enumerate(hs)
        hist_edges(h) == edges || throw(ArgumentError(
            "merge input $k has different bin edges; rebin the histograms to common " *
            "edges with rebin(h, edges) before merging"))
    end
    w, v = full_arrays(a)
    for h in hs[2:end]
        hw, hv = full_arrays(h)
        w .+= hw
        v .+= hv
    end
    all(isfinite, w) && all(isfinite, v) ||
        throw(ArgumentError("merged histogram contents overflow Float64"))
    return from_full(a, w, v, any(has_flow_arrays, hs), merged_entries(hs), merged_moments(hs))
end

Base.merge(a::AnyHist, b::AnyHist, more::AnyHist...) =
    throw(ArgumentError("merge needs histograms of the same type (Hist1D, Hist2D, or Hist3D)"))

# ---------------------------------------------------------------------------
# Projections

# Positions of each axis' (Σw·a, Σw·a²) in TH2/TH3::GetStats order, and the cross terms.
const MOMENT_AXIS = ((3, 4), (5, 6), (8, 9))
moment_cross(i::Int, j::Int) = (min(i, j), max(i, j)) == (1, 2) ? 7 :
                               (min(i, j), max(i, j)) == (1, 3) ? 10 : 11

project_moments(m, i::Int) =
    m === nothing ? nothing : (m[1], m[2], m[MOMENT_AXIS[i][1]], m[MOMENT_AXIS[i][2]])

project_moments(m, i::Int, j::Int) =
    m === nothing ? nothing :
    (m[1], m[2], m[MOMENT_AXIS[i][1]], m[MOMENT_AXIS[i][2]],
     m[MOMENT_AXIS[j][1]], m[MOMENT_AXIS[j][2]], m[moment_cross(i, j)])

"""True when a cell visible on the `kept` axes but in a flow bin of another axis is nonzero."""
function integrated_flow_hit(w, v, kept)
    shape = size(w)
    for index in CartesianIndices(w)
        c = Tuple(index)
        all(k -> 1 < c[k] < shape[k], kept) || continue
        any(k -> !(k in kept) && (c[k] == 1 || c[k] == shape[k]), eachindex(c)) || continue
        (w[index] != 0 || v[index] != 0) && return true
    end
    return false
end

"""
    projection(h::Union{Hist2D,Hist3D}, axis::Symbol; flow=false) -> Hist1D
    projection(h::Hist3D, (axis1, axis2); flow=false) -> Hist2D

Sum `Σw` and `Σw²` over the other axes. A 3D → 2D projection keeps the order
of the requested axes: `projection(h, (:z, :x))` has z on its x axis (ROOT
`Project3D("xz")`).

`flow=false` (default, unchanged): only visible bins are summed; the result has
no flow bins and unknown `entries` and moments, so ROOT shows bin-centre
statistics with effective entries. This differs from ROOT's default
`TH2::ProjectionX()`, which also sums the integrated axis' flow bins.

`flow=true`: integrated axes are summed **including** their flow bins (ROOT's
default projection range) and the projected axes keep their under/overflow
bins, so the result is a direct fill of the projected coordinates of all
events. Statistics are kept only when that is exact: if every cell that lands
in a visible result bin through an integrated-axis flow bin has `Σw == 0`
**and** `Σw² == 0`, the original `entries` and the matching unbinned moments
(with the cross term for 3D → 2D, reordered with the axes) are reused;
otherwise both are unknown (ROOT's fallback: bin-centre moments and effective
entries).

A `Hist2D`/`Hist3D` without `bins_with_flow` (built from visible arrays) is
taken to have **empty** flow bins: `flow=true` then sums zeros for them, the
result has zero under/overflow, and its stored `entries`/moments are reused.
This is exact only if the source really had no flow content; the flows cannot
be recovered from visible bins. Histograms from [`histogram`](@ref) store their
flow bins, so for them the result is exact.

Differences from ROOT 6.40.04 `ProjectionX`/`Project3D` (checked with plain
macros): over its default full range, ROOT reuses the statistics only when the
projection's total content, its flow bins included, equals the original
in-range `Σw` (relative 1e-12).

- Projected-axis flows: ROOT then falls back to bin-centre moments and
  effective entries; here the exact statistics are kept, equal to a direct
  fill of the projected coordinates.
- Signed weights cancelling in an integrated-axis flow cell (`Σw = 0`,
  `Σw² > 0`): ROOT keeps its statistics; here they are dropped, because that
  cell still adds `Σw²` (and events) to visible bins.
- In the fallback, ROOT's projection sets the in-range effective entries; the
  drawn result here gets `TH1::ResetStats` effective entries, which in ROOT
  6.40 include the result's flow bins. Mean and standard deviation agree.

See docs/recipes.md.
"""
function projection(h::Union{Hist2D,Hist3D}, axis::Symbol; flow::Bool=false)
    validate_hist_arrays(h)
    edges = hist_nd_edges(h)
    i = histogram_axis(axis, length(edges))
    dims = Tuple(j for j in 1:length(edges) if j != i)
    flow || return Hist1D(edges[i], vec(sum(h.sumw; dims)), vec(sum(h.sumw2; dims)))
    w, v = full_arrays(h)
    pw, pv = vec(sum(w; dims)), vec(sum(v; dims))
    all(isfinite, pw) && all(isfinite, pv) || throw(ArgumentError("projection overflows Float64"))
    clean = !integrated_flow_hit(w, v, (i,))
    return Hist1D(edges[i], pw[2:end-1], pv[2:end-1]; underflow=(pw[1], pv[1]),
                  overflow=(pw[end], pv[end]), entries=clean ? h.entries : nothing,
                  moments=clean ? project_moments(h.moments, i) : nothing)
end

function projection(h::Hist3D, axes::Tuple{Symbol,Symbol}; flow::Bool=false)
    validate_hist_arrays(h)
    i, j = map(a -> histogram_axis(a, 3), axes)
    i != j || throw(ArgumentError("projection axes must be distinct"))
    k = only(filter(a -> a != i && a != j, (1, 2, 3)))
    edges = hist_nd_edges(h)
    reduce2(a) = (r = dropdims(sum(a; dims=k); dims=k); i > j ? permutedims(r) : r)
    flow || return Hist2D(edges[i], edges[j], reduce2(h.sumw), reduce2(h.sumw2))
    w, v = full_arrays(h)
    pw, pv = reduce2(w), reduce2(v)
    all(isfinite, pw) && all(isfinite, pv) || throw(ArgumentError("projection overflows Float64"))
    clean = !integrated_flow_hit(w, v, (i, j))
    visible = (2:length(edges[i]), 2:length(edges[j]))
    return Hist2D(edges[i], edges[j], pw[visible...], pv[visible...];
                  bins_with_flow=has_flow_arrays(h) ? (sumw=pw, sumw2=pv) : nothing,
                  entries=clean ? h.entries : nothing,
                  moments=clean ? project_moments(h.moments, i, j) : nothing)
end
