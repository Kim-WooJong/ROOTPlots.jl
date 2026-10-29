function histogram_edges(edges::AbstractVector{<:Real})
    e = collect(Float64, edges)
    length(e) >= 2 && all(isfinite, e) && all(diff(e) .> 0) ||
        throw(ArgumentError("edges must be finite and increase strictly"))
    return e
end

function validate_hist_nd_arrays(edges::Tuple, sumw, sumw2)
    expected = map(e -> length(e) - 1, edges)
    all(e -> length(e) >= 2 && all(isfinite, e) && all(diff(e) .> 0), edges) &&
        size(sumw) == expected && size(sumw2) == expected &&
        all(isfinite, sumw) && all(isfinite, sumw2) && all(>=(0), sumw2) ||
        throw(ArgumentError("invalid histogram edges, bin dimensions, contents, or variances"))
    return expected
end

"""Owned visible-bin 2D histogram. Array dimensions are `(x, y)`; flow bins are not stored."""
struct Hist2D
    xedges::Vector{Float64}
    yedges::Vector{Float64}
    sumw::Matrix{Float64}
    sumw2::Matrix{Float64}
    function Hist2D(xedges::AbstractVector{<:Real}, yedges::AbstractVector{<:Real},
                    sumw::AbstractMatrix{<:Real}, sumw2::AbstractMatrix{<:Real})
        x, y = histogram_edges(xedges), histogram_edges(yedges)
        w, v = collect(Float64, sumw), collect(Float64, sumw2)
        validate_hist_nd_arrays((x, y), w, v)
        new(x, y, w, v)
    end
end

"""Owned visible-bin 3D histogram. Array dimensions are `(x, y, z)`; no flow bins are stored."""
struct Hist3D
    xedges::Vector{Float64}
    yedges::Vector{Float64}
    zedges::Vector{Float64}
    sumw::Array{Float64,3}
    sumw2::Array{Float64,3}
    function Hist3D(xedges::AbstractVector{<:Real}, yedges::AbstractVector{<:Real},
                    zedges::AbstractVector{<:Real}, sumw::AbstractArray{<:Real,3},
                    sumw2::AbstractArray{<:Real,3})
        x, y, z = histogram_edges(xedges), histogram_edges(yedges), histogram_edges(zedges)
        w, v = collect(Float64, sumw), collect(Float64, sumw2)
        validate_hist_nd_arrays((x, y, z), w, v)
        new(x, y, z, w, v)
    end
end

hist_nd_edges(h::Hist2D) = (h.xedges, h.yedges)
hist_nd_edges(h::Hist3D) = (h.xedges, h.yedges, h.zedges)

"""Revalidate owned but mutable histogram arrays before using their contents."""
validate_stats_hist_nd(h::Union{Hist2D,Hist3D}) =
    validate_hist_nd_arrays(hist_nd_edges(h), h.sumw, h.sumw2)

function validate_histogram_input(values::Tuple, weights, threads::Symbol)
    threads in (:auto, :on, :off) || throw(ArgumentError("threads must be :auto, :on, or :off"))
    n = length(first(values))
    all(v -> length(v) == n, values) || throw(ArgumentError("coordinate lengths must agree"))
    weights === nothing || weights isa AbstractVector ||
        throw(ArgumentError("weights must be a vector or nothing"))
    weights === nothing || length(weights) == n || throw(ArgumentError("weight length must match values"))
    Base.require_one_based_indexing(values...)
    weights === nothing || Base.require_one_based_indexing(weights)
    return n
end

"""Fill a vector using deterministic Julia CPU chunks and ROOT `[low, high)` bins.

Returns `(hist, underflow, overflow)`; each flow pair contains `(sumw, sumw2)`.
Signed weights are preserved, without normalization.
"""
function histogram(values::AbstractVector, edges::AbstractVector{<:Real};
                   weights=nothing, threads::Symbol=:auto)
    validate_histogram_input((values,), weights, threads)
    e = histogram_edges(edges)
    w, v = fill_chunks(values, weights, e, threads)
    return (; hist=Hist1D(e, w[2:end-1], v[2:end-1]),
            underflow=(w[1], v[1]), overflow=(w[end], v[end]))
end

function hist_nd_chunk_count(n::Int, shape::Tuple)
    # BigInt arithmetic prevents overflow in dimension products before allocation.
    bytes = 16 * prod(big.(shape))
    scratch_cap = Int(min(big(HIST_MAX_CHUNKS), max(big(1), fld(big(HIST_SCRATCH_BYTES), bytes))))
    return min(max(1, cld(n, HIST_CHUNK_SIZE)), HIST_MAX_CHUNKS, scratch_cap)
end

function fill_hist_nd_chunk!(sums, variances, values::Tuple, weights, edges::Tuple,
                             firstrow::Int, lastrow::Int)
    for row in firstrow:lastrow
        raw_x = map(v -> v[row], values)
        raw_w = event_weight(weights, row)
        all(x -> x isa Real && isfinite(x), raw_x) && raw_w isa Real && isfinite(raw_w) ||
            return (row, "non-finite or non-numeric event value/weight")
        x, w = map(Float64, raw_x), Float64(raw_w)
        all(isfinite, x) && isfinite(w) ||
            return (row, "event value/weight is not finite in Float64")
        variance = w * w
        isfinite(variance) || return (row, "event weight squared overflows Float64")
        slot = map((e, a) -> searchsortedlast(e, a) + 1, edges, x)
        new_sum, new_var = sums[slot...] + w, variances[slot...] + variance
        isfinite(new_sum) && isfinite(new_var) ||
            return (row, "histogram accumulator overflows Float64")
        sums[slot...], variances[slot...] = new_sum, new_var
    end
    return nothing
end

function fill_hist_nd_chunks(values::Tuple, weights, edges::Tuple, mode::Symbol)
    n = length(first(values))
    shape = map(e -> length(e) + 1, edges)
    nchunks = hist_nd_chunk_count(n, shape)
    sums = [zeros(shape) for _ in 1:nchunks]
    variances = [zeros(shape) for _ in 1:nchunks]
    problems = Vector{Union{Nothing,Tuple{Int,String}}}(nothing, nchunks)
    work(j) = (problems[j] = fill_hist_nd_chunk!(sums[j], variances[j], values, weights, edges,
                                                fld((j - 1) * n, nchunks) + 1, fld(j * n, nchunks)))
    parallel = mode != :off && Threads.nthreads(:default) > 1 && nchunks > 1 &&
               (mode == :on || n >= 100_000)
    if parallel
        Threads.@threads :dynamic for j in 1:nchunks
            work(j)
        end
    else
        for j in 1:nchunks
            work(j)
        end
    end
    for problem in problems
        problem === nothing || throw(ArgumentError("row $(problem[1]): $(problem[2])"))
    end
    # Reuse chunk one as the output, preserving a fixed chunk-major reduction.
    total, total_var = sums[1], variances[1]
    for j in 2:nchunks, k in eachindex(total)
        total[k] += sums[j][k]
        total_var[k] += variances[j][k]
        isfinite(total[k]) && isfinite(total_var[k]) ||
            throw(ArgumentError("merged histogram accumulator overflows Float64"))
    end
    return total, total_var
end

"""Fill 2D or 3D coordinates without ROOT objects or files.

`hist` contains visible bins. `bins_with_flow=(sumw, sumw2)` contains all bins,
in `(x, y[, z])` order: index 1 is underflow, indices `2:end-1` are visible,
and the final index is overflow on each axis, including mixed-axis flow bins.
The upper edge is overflow. `threads=:auto|:on|:off` changes scheduling only;
chunk boundaries and reduction order are deterministic and identical.
Private chunk storage is capped at 64 MiB, or one chunk when a single output
itself exceeds that size. No density normalization is applied.
"""
function histogram(values::Tuple{Vararg{AbstractVector,N}},
                   edges::Tuple{Vararg{AbstractVector,N}};
                   weights=nothing, threads::Symbol=:auto) where N
    N in (2, 3) || throw(ArgumentError("expected two or three coordinate vectors and edge vectors"))
    validate_histogram_input(values, weights, threads)
    all(e -> eltype(e) <: Real, edges) || throw(ArgumentError("edges must contain real numbers"))
    e = map(histogram_edges, edges)
    w, v = fill_hist_nd_chunks(values, weights, e, threads)
    visible = map(a -> 2:length(a), e)
    h = N == 2 ? Hist2D(e..., w[visible...], v[visible...]) :
                 Hist3D(e..., w[visible...], v[visible...])
    return (; hist=h, bins_with_flow=(sumw=w, sumw2=v))
end

function histogram_axis(axis::Symbol, n::Int)
    i = findfirst(==(axis), (:x, :y, :z))
    i === nothing || i > n ? throw(ArgumentError("invalid histogram projection axis $axis")) : i
end

"""Sum other visible dimensions and variances; flow bins never enter a projection."""
function projection(h::Union{Hist2D,Hist3D}, axis::Symbol)
    validate_stats_hist_nd(h)
    edges = hist_nd_edges(h)
    i = histogram_axis(axis, length(edges))
    dims = Tuple(j for j in 1:length(edges) if j != i)
    return Hist1D(edges[i], vec(sum(h.sumw; dims)), vec(sum(h.sumw2; dims)))
end

"""Project 3D visible bins onto the two requested axes, preserving their order."""
function projection(h::Hist3D, axes::Tuple{Symbol,Symbol})
    validate_stats_hist_nd(h)
    i, j = map(a -> histogram_axis(a, 3), axes)
    i != j || throw(ArgumentError("projection axes must be distinct"))
    k = only(filter(a -> a != i && a != j, (1, 2, 3)))
    w = dropdims(sum(h.sumw; dims=k); dims=k)
    v = dropdims(sum(h.sumw2; dims=k); dims=k)
    if i > j
        w, v = permutedims(w), permutedims(v)
    end
    edges = hist_nd_edges(h)
    return Hist2D(edges[i], edges[j], w, v)
end
