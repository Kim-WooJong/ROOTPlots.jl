# ROOT statistics and fit boxes (TPaveStats).
#
# The Julia side only normalizes options to ROOT's OptStat/OptFit integers and
# box placement; ROOT itself computes and prints every value.

"""Statistics-box fields: `(decimal digit position, digit value)` of ROOT's `OptStat`."""
const STATS_FIELDS = (name=(1,1), entries=(2,1), mean=(3,1), mean_error=(3,2),
                      std=(4,1), std_error=(4,2), underflow=(5,1), overflow=(6,1),
                      integral=(7,1), integral_width=(7,2), skewness=(8,1),
                      skewness_error=(8,2), kurtosis=(9,1), kurtosis_error=(9,2))
const STATS_LETTERS = Dict('n'=>(1,1), 'e'=>(2,1), 'm'=>(3,1), 'M'=>(3,2), 'r'=>(4,1),
                           'R'=>(4,2), 'u'=>(5,1), 'o'=>(6,1), 'i'=>(7,1), 'I'=>(7,2),
                           's'=>(8,1), 'S'=>(8,2), 'k'=>(9,1), 'K'=>(9,2))
"""Fit-box fields: `(decimal digit position, digit value)` of ROOT's `OptFit`."""
const FIT_FIELDS = (values=(1,1), free_values=(1,2), errors=(2,1), chi2=(3,1), prob=(4,1))
const FIT_LETTERS = Dict('v'=>(1,1), 'e'=>(2,1), 'c'=>(3,1), 'p'=>(4,1))
const STATS_POSITIONS = (:auto, :topright, :topleft, :bottomright, :bottomleft)

function stats_digits_code(pairs, ndigits)
    code = zeros(Int, ndigits)
    for (position, value) in pairs
        code[position] = max(code[position], value)
    end
    return sum(code[i] * 10^(i-1) for i in 1:ndigits)
end

function stats_code(value, fields, letters, ndigits, default, one, what)
    (value === false || value === nothing) && return 0
    value === true && return default
    if value isa Integer
        0 <= value < 10^ndigits && all(<=(2), digits(value)) ||
            throw(ArgumentError("integer $what must use at most $ndigits ROOT digits, each 0, 1, or 2"))
        return value == 1 ? one : Int(value) # ROOT reads 1 as its default content
    elseif value isa AbstractString
        isempty(value) && return 0
        all(c -> haskey(letters, c), value) ||
            throw(ArgumentError("ROOT $what letters must be from \"$(join(sort(collect(keys(letters)))))\""))
        return stats_digits_code((letters[c] for c in value), ndigits)
    end
    items = value isa Symbol ? (value,) : value
    items isa Union{Tuple,AbstractVector} && !isempty(items) && all(x -> x isa Symbol, items) ||
        throw(ArgumentError("$what must be Bool, a ROOT integer/letter string, or Symbols"))
    for item in items
        haskey(fields, item) ||
            throw(ArgumentError("unknown $what field :$item; valid fields: $(join(keys(fields), ", "))"))
    end
    return stats_digits_code((fields[item] for item in items), ndigits)
end

"""
    stats_option(stats) -> Int

Normalize a `stats` keyword to ROOT's `gStyle->SetOptStat` integer.
`true` is ROOT's default `1111` (name, entries, mean, std dev); `false` is 0.
Also accepted: a tuple of fields (`:name, :entries, :mean, :mean_error, :std,
:std_error, :underflow, :overflow, :integral, :integral_width, :skewness,
:skewness_error, :kurtosis, :kurtosis_error`), a ROOT letter string such as
`"ksiourmen"` (upper case adds errors), or a ROOT integer such as `1001111`
(the integer 1 means ROOT's default 1111).
"""
stats_option(stats) = stats_code(stats, STATS_FIELDS, STATS_LETTERS, 9, 1111, 1111, "stats")

"""
    fit_option(fitstats) -> Int

Normalize a `fitstats` keyword to ROOT's `gStyle->SetOptFit` integer. `true` is
`1111` (probability, χ²/ndf, errors, values). Fields are `:prob, :chi2, :errors,
:values, :free_values`; ROOT letters `"pcev"` and integers are also accepted
(the integer 1 means ROOT's `111`).
"""
fit_option(fitstats) = stats_code(fitstats, FIT_FIELDS, FIT_LETTERS, 4, 1111, 111, "fitstats")

"""Normalized statistics/fit box request for one ROOT pad (internal)."""
struct StatsRequest
    optstat::Int
    optfit::Int
    corner::String
    box::NTuple{4,Float64}
    name::String
    headroom::Bool
end

stats_active(r::StatsRequest) = r.optstat > 0 || r.optfit > 0

function stats_box_position(position)
    if position isa Symbol
        position in STATS_POSITIONS ||
            throw(ArgumentError("stats_position must be one of $(STATS_POSITIONS) or NDC (x1, y1, x2, y2)"))
        return String(position), (0.0, 0.0, 0.0, 0.0)
    end
    position isa Union{Tuple,AbstractVector} && length(position) == 4 &&
        all(x -> x isa Real && isfinite(x), position) ||
        throw(ArgumentError("stats_position must be a corner Symbol or NDC (x1, y1, x2, y2)"))
    x1, y1, x2, y2 = Float64.(Tuple(position))
    0 <= x1 < x2 <= 1 && 0 <= y1 < y2 <= 1 ||
        throw(ArgumentError("stats_position NDC box needs 0 ≤ x1 < x2 ≤ 1 and 0 ≤ y1 < y2 ≤ 1"))
    return "manual", (x1, y1, x2, y2)
end

"""
Build the box request. `name` is the header line (ROOT prints the object
name); `fit` is the attached `FitResult` or `nothing`; `allow_stats=false`
marks objects without ROOT histogram statistics (graphs, stacks).
`stats_headroom=true` asks the painter to raise the automatic y maximum of a
1-D frame so the data stay below the boxes; it conflicts with explicit
`ylimits` and needs a 1-D frame (`one_dimensional`).
"""
function stats_request(; stats=false, fitstats=false, stats_position=:auto, name="Histogram",
                       fit=nothing, allow_stats=true, what="this plot",
                       stats_headroom=false, ylimits=nothing, one_dimensional=true)
    optstat = stats_option(stats)
    optfit = fit_option(fitstats)
    optstat > 0 && !allow_stats &&
        throw(ArgumentError("ROOT draws no statistics box for $what"))
    optfit > 0 && fit === nothing &&
        throw(ArgumentError("fitstats needs attached_fit=FitResult (the fit shown in the box)"))
    corner, box = stats_box_position(stats_position)
    stats_headroom isa Bool || throw(ArgumentError("stats_headroom must be true or false"))
    stats_headroom && !one_dimensional &&
        throw(ArgumentError("stats_headroom applies to 1-D frames (histograms, profiles, graphs)"))
    stats_headroom && ylimits !== nothing &&
        throw(ArgumentError("stats_headroom raises the automatic y maximum; it cannot be combined with ylimits"))
    label = String(name)
    clean_field(label)
    occursin('\0', label) && throw(ArgumentError("NUL is not permitted in stats_name"))
    return StatsRequest(optstat, optfit, corner, box, isempty(label) ? "Histogram" : label, stats_headroom)
end

function write_stats_row(io, r::StatsRequest)
    stats_active(r) || return nothing
    println(io, join(("STATSBOX", r.optstat, r.optfit, r.corner, r.box...,
                      clean_field(root_latex(r.name)), Int(r.headroom)), '\t'))
    return nothing
end

"""Reject fits that cannot be attached to a native TF1 (status, model, parameters)."""
function validate_attached_fit(fit)
    fit isa FitResult && fit.status == 0 && fit.valid && haskey(ROOT_FIT_NPAR, fit.model) &&
        length(fit.parameters) == ROOT_FIT_NPAR[fit.model] == length(fit.errors) == length(fit.fixed) &&
        all(isfinite, fit.parameters) && all(e -> isfinite(e) && e >= 0, fit.errors) &&
        isfinite(fit.chi2) && fit.chi2 >= 0 && fit.ndf >= 0 &&
        all(isfinite, fit.range) && fit.range[1] < fit.range[2] ||
        throw(ArgumentError("attached_fit must be a valid supported ROOT fit"))
    validate_model_parameters(fit.model, fit.parameters)
    return fit
end

"""Serialize a fitted TF1: parameters, then errors, χ², ndf, and fixed flags for the fit box."""
function write_attached_fit(io, fit::FitResult)
    println(io, join(("FITFUNC", String(fit.model), fit.range..., length(fit.parameters),
                      fit.parameters...), '\t'))
    println(io, join(("FITERR", fit.chi2, fit.ndf, fit.errors..., Int.(fit.fixed)...), '\t'))
end

# Histogram bookkeeping for ROOT: entries, flow bins (ROOT indices), unbinned moments.
function write_hist_stats(io, h::Hist1D)
    n = length(h.sumw)
    m = h.moments === nothing ? () : h.moments
    println(io, join(("HSTATS", h.entries === nothing ? "-" : h.entries, 2, length(m), m...), '\t'))
    println(io, join(("FLOW", 0, h.underflow...), '\t'))
    println(io, join(("FLOW", n + 1, h.overflow...), '\t'))
end

function write_hist_stats(io, h::Union{Hist2D,Hist3D})
    m = h.moments === nothing ? () : h.moments
    cells = Tuple{Tuple,Float64,Float64}[]
    if h.bins_with_flow !== nothing
        w, v = h.bins_with_flow.sumw, h.bins_with_flow.sumw2
        shape = size(w)
        for index in CartesianIndices(w)
            i = Tuple(index)
            any(k -> i[k] == 1 || i[k] == shape[k], eachindex(i)) || continue
            (w[index] != 0 || v[index] != 0) && push!(cells, (i .- 1, w[index], v[index]))
        end
    end
    println(io, join(("HSTATS", h.entries === nothing ? "-" : h.entries, length(cells),
                      length(m), m...), '\t'))
    for (i, w, v) in cells
        println(io, join(("FLOW", i..., w, v), '\t'))
    end
end

"""Pick the box header: explicit `stats_name`, else the label, else a type name."""
stats_header(stats_name, label, fallback="Histogram") =
    stats_name !== nothing ? String(stats_name) :
    !isempty(String(label)) ? String(label) : fallback
