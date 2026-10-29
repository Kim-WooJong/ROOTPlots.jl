"""
    Hist1D(edges, sumw, sumw2; underflow=(0, 0), overflow=(0, 0), entries=nothing, moments=nothing)

Visible-bin histogram with ROOT-style half-open edges. Inputs are edges, Σw and
Σw² of the visible bins; no file or ROOT object is owned.

The optional bookkeeping mirrors what a ROOT `TH1D` filled with `Fill` knows and
is used only when the histogram is drawn by native ROOT (e.g. in a statistics box):

- `underflow`, `overflow`: `(Σw, Σw²)` of ROOT bins 0 and n+1.
- `entries`: the number of fill calls (ROOT `GetEntries()`), including flow
  entries, or `nothing` when unknown.
- `moments`: ROOT's unbinned statistics `(Σw, Σw², Σwx, Σwx²)` of the
  in-range entries (`TH1::GetStats` order), or `nothing`.

[`histogram`](@ref) sets all of them. A histogram built directly from bin
contents keeps the defaults: zero flows, unknown entries (ROOT then shows the
effective entries `(Σw)²/Σw²`, as `TH1::ResetStats` does), and bin-centre moments.
"""
struct Hist1D
    edges::Vector{Float64}
    sumw::Vector{Float64}
    sumw2::Vector{Float64}
    underflow::Tuple{Float64,Float64}
    overflow::Tuple{Float64,Float64}
    entries::Union{Nothing,Int}
    moments::Union{Nothing,NTuple{4,Float64}}

    function Hist1D(edges::AbstractVector{<:Real}, sumw::AbstractVector{<:Real},
                    sumw2::AbstractVector{<:Real}; underflow=(0.0, 0.0), overflow=(0.0, 0.0),
                    entries=nothing, moments=nothing)
        e, w, v = Float64.(edges), Float64.(sumw), Float64.(sumw2)
        length(e) >= 2 && length(w) == length(e) - 1 && length(v) == length(w) ||
            throw(ArgumentError("expected n+1 edges and n visible bin values/variances"))
        all(isfinite, e) && all(diff(e) .> 0) || throw(ArgumentError("edges must increase strictly"))
        all(isfinite, w) && all(isfinite, v) && all(v .>= 0) ||
            throw(ArgumentError("non-finite bin or invalid variance"))
        new(e, w, v, flow_pair(underflow), flow_pair(overflow), entry_count(entries),
            hist_moments(moments, 4))
    end
end

"""Validate one `(Σw, Σw²)` flow pair."""
function flow_pair(flow)
    flow isa Tuple && length(flow) == 2 && all(x -> x isa Real && isfinite(x), flow) &&
        flow[2] >= 0 || throw(ArgumentError("flow must be a finite (Σw, Σw²) pair with Σw² ≥ 0"))
    return (Float64(flow[1]), Float64(flow[2]))
end

"""Validate an optional fill count (ROOT `GetEntries`)."""
function entry_count(entries)
    entries === nothing && return nothing
    entries isa Integer && entries >= 0 ||
        throw(ArgumentError("entries must be a nonnegative integer or nothing"))
    return Int(entries)
end

"""Validate optional unbinned ROOT statistics (`TH1::GetStats` order)."""
function hist_moments(moments, n::Int)
    moments === nothing && return nothing
    length(moments) == n && all(x -> x isa Real && isfinite(x), moments) ||
        throw(ArgumentError("moments must hold $n finite values in ROOT TH1::GetStats order"))
    m = ntuple(i -> Float64(moments[i]), n)
    m[2] >= 0 || throw(ArgumentError("moment Σw² must be nonnegative"))
    return m
end

const HIST_CHUNK_SIZE = 65_536
const HIST_MAX_CHUNKS = 64
const HIST_SCRATCH_BYTES = 64 * 1024 * 1024

event_weight(::Nothing, ::Int) = 1.0
event_weight(col::AbstractVector, i::Int) = col[i]

"""Fill one disjoint event range, returning its first invalid row without a worker TaskFailedException.

`moments` accumulates ROOT's in-range statistics `(Σw, Σw², Σwx, Σwx²)`.
"""
function fill_chunk!(sums::Vector{Float64}, variances::Vector{Float64},
                     moments::Vector{Float64},
                     values::AbstractVector, weights, edges::Vector{Float64},
                     firstrow::Int, lastrow::Int)
    nslots = length(sums)
    for row in firstrow:lastrow
        raw_x = values[row]
        raw_w = event_weight(weights, row)
        if !(raw_x isa Real && raw_w isa Real && isfinite(raw_x) && isfinite(raw_w))
            return (row, "non-finite or non-numeric event value/weight")
        end
        x, w = Float64(raw_x), Float64(raw_w)
        if !(isfinite(x) && isfinite(w))
            return (row, "event value/weight is not finite in Float64")
        end
        variance = w * w
        isfinite(variance) || return (row, "event weight squared overflows Float64")
        slot = searchsortedlast(edges, x) + 1
        new_sum = sums[slot] + w
        new_var = variances[slot] + variance
        if !(isfinite(new_sum) && isfinite(new_var))
            return (row, "histogram accumulator overflows Float64")
        end
        sums[slot], variances[slot] = new_sum, new_var
        if 1 < slot < nslots
            # Like TH1::Fill, only in-range entries enter the unbinned moments.
            moments[1] += w; moments[2] += variance
            moments[3] += w * x; moments[4] += w * x * x
        end
    end
    return nothing
end

function fill_chunks(values::AbstractVector, weights, edges::Vector{Float64}, mode::Symbol)
    n = length(values)
    nslots = length(edges) + 1 # underflow, visible bins, overflow
    scratch_cap = max(1, fld(HIST_SCRATCH_BYTES, 16 * nslots))
    nchunks = min(max(1, cld(n, HIST_CHUNK_SIZE)), HIST_MAX_CHUNKS, scratch_cap)
    sums = [zeros(nslots) for _ in 1:nchunks]
    variances = [zeros(nslots) for _ in 1:nchunks]
    moments = [zeros(4) for _ in 1:nchunks]
    problems = Vector{Union{Nothing,Tuple{Int,String}}}(nothing, nchunks)
    work(j) = (problems[j] = fill_chunk!(sums[j], variances[j], moments[j], values, weights, edges,
                                        fld((j - 1) * n, nchunks) + 1,
                                        fld(j * n, nchunks)))
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
    total = zeros(nslots)
    total_var = zeros(nslots)
    for j in 1:nchunks, k in 1:nslots
        total[k] += sums[j][k]
        total_var[k] += variances[j][k]
        isfinite(total[k]) && isfinite(total_var[k]) ||
            throw(ArgumentError("merged histogram accumulator overflows Float64"))
    end
    return total, total_var, merge_moments(moments)
end

"""Reduce per-chunk moments in chunk order; `nothing` if they overflow Float64."""
function merge_moments(chunks)
    total = zeros(length(first(chunks)))
    for m in chunks
        total .+= m
    end
    return all(isfinite, total) ? Tuple(total) : nothing
end
