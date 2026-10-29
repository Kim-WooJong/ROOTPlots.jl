"""Visible-bin histogram with ROOT-style half-open edges. Inputs are edges, Σw and Σw²; no flow bins, file, or ROOT object are owned."""
struct Hist1D
    edges::Vector{Float64}
    sumw::Vector{Float64}
    sumw2::Vector{Float64}

    function Hist1D(edges::AbstractVector{<:Real}, sumw::AbstractVector{<:Real},
                    sumw2::AbstractVector{<:Real})
        e, w, v = Float64.(edges), Float64.(sumw), Float64.(sumw2)
        length(e) >= 2 && length(w) == length(e) - 1 && length(v) == length(w) ||
            throw(ArgumentError("expected n+1 edges and n visible bin values/variances"))
        all(isfinite, e) && all(diff(e) .> 0) || throw(ArgumentError("edges must increase strictly"))
        all(isfinite, w) && all(isfinite, v) && all(v .>= 0) ||
            throw(ArgumentError("non-finite bin or invalid variance"))
        new(e, w, v)
    end
end

const HIST_CHUNK_SIZE = 65_536
const HIST_MAX_CHUNKS = 64
const HIST_SCRATCH_BYTES = 64 * 1024 * 1024

event_weight(::Nothing, ::Int) = 1.0
event_weight(col::AbstractVector, i::Int) = col[i]

"""Fill one disjoint event range, returning its first invalid row without a worker TaskFailedException."""
function fill_chunk!(sums::Vector{Float64}, variances::Vector{Float64},
                     values::AbstractVector, weights, edges::Vector{Float64},
                     firstrow::Int, lastrow::Int)
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
    problems = Vector{Union{Nothing,Tuple{Int,String}}}(nothing, nchunks)
    work(j) = (problems[j] = fill_chunk!(sums[j], variances[j], values, weights, edges,
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
    return total, total_var
end
