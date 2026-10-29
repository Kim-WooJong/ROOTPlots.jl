"""Native ROOT graph line/marker/fill style. Fill defaults to the line color at 35% opacity."""
Base.@kwdef struct GraphStyle
    color::String = "#000000"
    marker_style::Int = 20
    marker_size::Float64 = 0.8
    line_width::Int = 2
    line_style::Int = 1
    fill_color::Union{Nothing,String} = nothing
    fill_style::Int = 1001
    fill_alpha::Float64 = 0.35
end

# The original five-field constructor remains valid.
GraphStyle(color, marker_style, marker_size, line_width, line_style) =
    GraphStyle(;color,marker_style,marker_size,line_width,line_style)

function check_graph_style(s::GraphStyle)
    occursin(r"^#[0-9a-fA-F]{6}$", s.color) || throw(ArgumentError("graph color must be #RRGGBB"))
    1 <= s.marker_style <= 50 || throw(ArgumentError("marker style must be 1:50"))
    isfinite(s.marker_size) && s.marker_size > 0 || throw(ArgumentError("marker size must be positive"))
    s.line_width > 0 && 1 <= s.line_style <= 10 || throw(ArgumentError("invalid line style/width"))
    (s.fill_color === nothing || occursin(r"^#[0-9a-fA-F]{6}$",s.fill_color)) ||
        throw(ArgumentError("graph fill color must be #RRGGBB or nothing"))
    0 <= s.fill_style <= 9999 && isfinite(s.fill_alpha) && 0 <= s.fill_alpha <= 1 ||
        throw(ArgumentError("invalid graph fill style/alpha"))
    return s
end

function graph_errors(input, n, axis)
    input === nothing && return (zeros(n), zeros(n), :plain)
    if input isa AbstractVector
        e = Float64.(input)
        kind = :symmetric
        lo, hi = e, copy(e)
    elseif input isa BinErrors || (input isa NamedTuple && keys(input) == (:low, :high))
        lo, hi = Float64.(input.low), Float64.(input.high)
        kind = :asymmetric
    else
        throw(ArgumentError("$axis errors need a vector, (low=…, high=…), or BinErrors"))
    end
    length(lo) == n == length(hi) && all(isfinite, lo) && all(isfinite, hi) &&
        all(>=(0), lo) && all(>=(0), hi) ||
        throw(ArgumentError("invalid $axis errors"))
    return (lo, hi, kind)
end

"""Owned XY graph points with absolute symmetric/asymmetric errors.
Connected lines retain input order. `draw` presets are `:points`, `:line`,
`:linepoints`, `:curve`, `:curvepoints`, `:fill`, `:bar`, and `:bar_offset`.
`errors` supports `:bars`, `:none`, `:band`, `:smoothband`, `:boxes`, `:caps`,
`:nocaps`, `:arrows`, and `:filledarrows`. Smoothed drawings are ROOT display
interpolation, not a numerical fit or additional uncertainty calculation.
A raw ROOT `draw` string overrides `errors`; the owned frame rejects A/SAME
and adds SAME itself, unlike histogram raw-option forwarding.
"""
struct GraphSeries
    x::Vector{Float64}
    y::Vector{Float64}
    exlow::Vector{Float64}
    exhigh::Vector{Float64}
    eylow::Vector{Float64}
    eyhigh::Vector{Float64}
    kind::Symbol
    label::String
    draw::Union{Symbol,String}
    errors::Symbol
    style::GraphStyle
end

const GRAPH_DRAW_MODES = (:points,:line,:linepoints,:curve,:curvepoints,:fill,:bar,:bar_offset)
const GRAPH_ERROR_MODES = (:bars,:none,:band,:smoothband,:boxes,:caps,:nocaps,:arrows,:filledarrows)

"""Validate per-layer ROOT graph options. The canvas owns axes and appends SAME.
Raw options override `errors`; put e.g. `3`, `Z`, or `X` in the option itself.
"""
function graph_root_option(option::AbstractString)
    value = root_draw_option(option)
    occursin('A',uppercase(value)) &&
        throw(ArgumentError("graph axes belong to GraphSpec; raw draw must not contain A or SAME"))
    return value
end

graph_visible_errors(s::GraphSeries) = s.kind != :plain &&
    (s.draw isa String ? !occursin('X',uppercase(s.draw)) : s.errors != :none)
graph_series_draw(s::GraphSeries) = s.draw isa String ? "ROOT:" * s.draw : String(s.draw)

function GraphSeries(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                     xerrors=nothing, yerrors=nothing, label="", draw=:points,
                     errors=:bars, style=GraphStyle())
    xx, yy = Float64.(x), Float64.(y)
    length(xx) == length(yy) && !isempty(xx) && all(isfinite, xx) && all(isfinite, yy) ||
        throw(ArgumentError("graph needs equal, nonempty, finite x/y arrays"))
    exl, exh, xkind = graph_errors(xerrors, length(xx), :x)
    eyl, eyh, ykind = graph_errors(yerrors, length(xx), :y)
    kind = :asymmetric in (xkind, ykind) ? :asymmetric :
           :symmetric in (xkind, ykind) ? :symmetric : :plain
    s = GraphSeries(xx, yy, exl, exh, eyl, eyh, kind, String(label),
                    draw isa AbstractString ? String(draw) : draw, errors, style)
    validate_graph(s)
    return s
end

function validate_graph(s::GraphSeries)
    n = length(s.x)
    n > 0 && all(length(v) == n for v in (s.y,s.exlow,s.exhigh,s.eylow,s.eyhigh)) ||
        throw(ArgumentError("graph arrays changed length"))
    all(isfinite, s.x) && all(isfinite, s.y) ||
        throw(ArgumentError("graph coordinates must be finite"))
    for e in (s.exlow,s.exhigh,s.eylow,s.eyhigh)
        all(isfinite, e) && all(>=(0), e) || throw(ArgumentError("invalid graph errors"))
    end
    all(isfinite, s.x .- s.exlow) && all(isfinite, s.x .+ s.exhigh) &&
        all(isfinite, s.y .- s.eylow) && all(isfinite, s.y .+ s.eyhigh) ||
        throw(ArgumentError("graph error endpoints overflow"))
    s.kind in (:plain,:symmetric,:asymmetric) || throw(ArgumentError("invalid graph kind"))
    if s.draw isa String
        graph_root_option(s.draw)
    else
        s.draw in GRAPH_DRAW_MODES || throw(ArgumentError("invalid graph draw mode"))
    end
    s.errors in GRAPH_ERROR_MODES || throw(ArgumentError("invalid graph error mode"))
    s.errors in (:band,:smoothband,:boxes,:caps,:nocaps,:arrows,:filledarrows) && s.kind == :plain &&
        throw(ArgumentError("the requested error display requires stored graph errors"))
    clean_field(s.label)
    occursin('\0', s.label) && throw(ArgumentError("NUL is not permitted in a graph label"))
    check_graph_style(s.style)
    return s
end

"""Single ROOT graph canvas with physical-coordinate axis limits."""
Base.@kwdef struct GraphSpec
    xlabel::String = "x"
    ylabel::String = "y"
    title::String = ""
    xlimits::Union{Nothing,Tuple{Float64,Float64}} = nothing
    ylimits::Union{Nothing,Tuple{Float64,Float64}} = nothing
    logx::Bool = false
    logy::Bool = false
    legend::Symbol = :topright
    size::Tuple{Int,Int} = (800,600)
    xaxis::AxisStyle = AxisStyle()
    yaxis::AxisStyle = AxisStyle()
    canvas::CanvasStyle = CanvasStyle()
    lines::Vector{PlotLine} = PlotLine[]
    texts::Vector{PlotText} = PlotText[]
    legend_spec::Union{Nothing,LegendSpec} = nothing
end

GraphSpec(xlabel,ylabel,title,xlimits,ylimits,logx,logy,legend,size) =
    GraphSpec(;xlabel,ylabel,title,xlimits,ylimits,logx,logy,legend,size)
with_graph_spec(spec::GraphSpec; kwargs...) = GraphSpec(;
    merge(NamedTuple{fieldnames(GraphSpec)}(Tuple(getfield(spec,f) for f in fieldnames(GraphSpec))),
          (;kwargs...))...)

function graph_limits(series::AbstractVector{GraphSeries}, spec::GraphSpec)
    validate_scene(spec; nseries=length(series))
    isempty(series) && throw(ArgumentError("at least one graph series is needed"))
    spec.legend in (:topright,:topleft,:bottomright,:bottomleft,:none) ||
        throw(ArgumentError("invalid legend position"))
    all(x -> 1 <= x <= 10000, spec.size) || throw(ArgumentError("graph size must be in 1:10000"))
    for label in (spec.xlabel,spec.ylabel,spec.title)
        clean_field(label)
        occursin('\0', label) && throw(ArgumentError("NUL is not permitted in labels"))
    end
    foreach(validate_graph, series)
    function axis(which, explicit, logarithmic)
        if logarithmic
            for s in series
                center = which == :x ? s.x : s.y
                elo = which == :x ? s.exlow : s.eylow
                minimum(center) > 0 || throw(ArgumentError("log $which requires positive graph centers"))
                !graph_visible_errors(s) || minimum(center .- elo) > 0 ||
                    throw(ArgumentError("log $which requires positive visible error endpoints"))
            end
        end
        if explicit !== nothing
            lo, hi = explicit
            isfinite(lo) && isfinite(hi) && lo < hi && (!logarithmic || lo > 0) ||
                throw(ArgumentError("invalid $which limits"))
            return (lo,hi)
        end
        lows = Float64[]
        highs = Float64[]
        for s in series
            center = which == :x ? s.x : s.y
            elo = which == :x ? s.exlow : s.eylow
            ehi = which == :x ? s.exhigh : s.eyhigh
            use_errors = graph_visible_errors(s)
            append!(lows, use_errors ? center .- elo : center)
            append!(highs, use_errors ? center .+ ehi : center)
            which == :y && !logarithmic && s.draw in (:bar,:bar_offset) && push!(lows,0.)
        end
        lo, hi = minimum(lows), maximum(highs)
        if logarithmic
            lo > 0 || throw(ArgumentError("log $which requires positive centers and visible error endpoints"))
            lo == hi && return (lo/sqrt(10), hi*sqrt(10))
            pad = (log(hi)-log(lo))*0.05
            return (exp(log(lo)-pad), exp(log(hi)+pad))
        end
        lo == hi && return (lo-max(0.05abs(lo),0.5), hi+max(0.05abs(hi),0.5))
        pad = 0.05*(hi-lo)
        return (lo-pad,hi+pad)
    end
    xlim = axis(:x,spec.xlimits,spec.logx)
    ylim = axis(:y,spec.ylimits,spec.logy)
    for (lim, logarithmic) in ((xlim,spec.logx),(ylim,spec.logy))
        all(isfinite,lim) && lim[1] < lim[2] && (!logarithmic || lim[1] > 0) ||
            throw(ArgumentError("automatic graph axis limits are not representable"))
    end
    return (xlim, ylim)
end

const ROOT_FIT_MODELS = Dict(:pol0=>"pol0", :pol1=>"pol1", :pol2=>"pol2",
    :gaus=>"gaus", :exp=>"[0]*exp([1]*x)",
    :erfc_edge=>"[0]+0.5*[1]*TMath::Erfc((x-[2])/(TMath::Sqrt(2.0)*[3]))",
    :double_erf=>"[0]+0.5*[1]*(TMath::Erf((x-[2])/(TMath::Sqrt(2.0)*[4]))-TMath::Erf((x-[3])/(TMath::Sqrt(2.0)*[5])))")
const ROOT_FIT_NPAR = Dict(:pol0=>1,:pol1=>2,:pol2=>3,:gaus=>3,:exp=>2,
                           :erfc_edge=>4,:double_erf=>6)
edge_width_indices(model) = model == :erfc_edge ? (4,) : model == :double_erf ? (5,6) : ()

function validate_model_parameters(model, p)
    haskey(ROOT_FIT_NPAR,model) && length(p) == ROOT_FIT_NPAR[model] && all(isfinite,p) ||
        throw(ArgumentError("invalid model parameters"))
    all(i -> p[i] > 0, edge_width_indices(model)) ||
        throw(ArgumentError("edge widths must be strictly positive"))
    return p
end

"""Numerical ROOT fit result; status/valid/covariance_status are ROOT's, not an independent quality guarantee."""
struct FitResult
    model::Symbol
    range::Tuple{Float64,Float64}
    parameters::Vector{Float64}
    errors::Vector{Float64}
    covariance::Matrix{Float64}
    fixed::BitVector
    status::Int
    valid::Bool
    covariance_status::Int
    chi2::Float64
    chi2_kind::Symbol
    ndf::Int
    npoints::Int
    edm::Float64
    minimizer::String
    normalized_errors::Bool
end

function model_value(model::Symbol, p, x)
    validate_model_parameters(model,p)
    model == :pol0 && return p[1]
    model == :pol1 && return p[1]+p[2]*x
    model == :pol2 && return p[1]+p[2]*x+p[3]*x*x
    model == :gaus && return p[1]*exp(-0.5*((x-p[2])/p[3])^2)
    model == :exp && return p[1]*exp(p[2]*x)
    model == :erfc_edge && return p[1]+0.5p[2]*erfc((x-p[3])/(sqrt(2)*p[4]))
    model == :double_erf && return p[1]+0.5p[2]*(erf((x-p[3])/(sqrt(2)*p[5]))-
                                               erf((x-p[4])/(sqrt(2)*p[6])))
    throw(ArgumentError("unsupported fit model"))
end

"""Sample a ROOT fit's returned parameters for a graph overlay; sampling is explicit, not a confidence band."""
function fit_curve(result::FitResult; n::Int=300, label="ROOT fit",
                   style=GraphStyle(color="#cc3311"), xrange=result.range)
    n >= 2 || throw(ArgumentError("fit curve needs at least two samples"))
    lo, hi = Float64.(xrange)
    isfinite(lo) && isfinite(hi) && lo < hi || throw(ArgumentError("invalid curve range"))
    xx = collect(range(lo, hi; length=n))
    yy = [model_value(result.model, result.parameters, x) for x in xx]
    return GraphSeries(xx, yy; label, style, draw=:line, errors=:none)
end

function graphfit_executable(config::AbstractString; force::Bool=false)
    source = normpath(joinpath(@__DIR__, "..", "backend", "root", "graphfit.cpp"))
    return compile_root_backend(config, source, "rootplots-graphfit", "-std=c++20";
                                force, dependencies=[joinpath(dirname(source), "advanced_histograms.hpp"),
                                                     joinpath(dirname(source), "layouts.hpp"),
                                                     joinpath(dirname(source), "primitives.hpp"),
                                                     joinpath(dirname(source), "specialgraphs.hpp")])
end

function write_graph_series(io, series)
    for (id,s) in enumerate(series)
        row = ("SERIES",id,length(s.x),s.kind,graph_series_draw(s),s.errors,
                          clean_field(root_latex(s.label)),s.style.color,
                          s.style.marker_style,s.style.marker_size,
                          s.style.line_width,s.style.line_style)
        if s.style.fill_color !== nothing || s.style.fill_style != 1001 || s.style.fill_alpha != 0.35
            row = (row...,something(s.style.fill_color,"-"),s.style.fill_style,s.style.fill_alpha)
        end
        println(io, join(row, '\t'))
        for i in eachindex(s.x)
            println(io, join(("POINT",id,i,s.x[i],s.y[i],s.exlow[i],s.exhigh[i],
                              s.eylow[i],s.eyhigh[i]), '\t'))
        end
    end
    println(io,"END\t",length(series))
end

function write_graph_spec(path, series, spec, limits)
    open(path, "w") do io
        println(io, "ROOTPlotsGraph\t1")
        println(io, join(("DRAW",spec.size[1],spec.size[2],Int(spec.logx),Int(spec.logy),
                          limits[1]...,limits[2]...,String(spec.legend)), '\t'))
        println(io, join(("LABELS",clean_field(root_latex(spec.title)),
                          clean_field(root_latex(spec.xlabel)),clean_field(root_latex(spec.ylabel))), '\t'))
        write_scene(io,spec)
        write_graph_series(io,series)
    end
end

"""Write one native two-pad χ² histogram-fit diagnostic request."""
function write_fit_plot_spec(path, series, spec, limits, pull_limits)
    open(path, "w") do io
        println(io, "ROOTPlotsGraph\t1")
        println(io, join(("DRAWFIT",spec.size[1],spec.size[2],Int(spec.logx),Int(spec.logy),
                          limits[1]...,limits[2]...,pull_limits...,String(spec.legend)), '\t'))
        println(io, join(("LABELS",clean_field(root_latex(spec.title)),
                          clean_field(root_latex(spec.xlabel)),clean_field(root_latex(spec.ylabel))), '\t'))
        write_graph_series(io,series)
    end
end

"""Render one-pad ROOT TGraph/TGraphErrors/TGraphAsymmErrors series in the selected format."""
function render_root(series::AbstractVector{GraphSeries}, spec::GraphSpec,
                     outstem::AbstractString; format=:png, root_config=nothing)
    paths = plot_output_paths(outstem, format)
    limits = graph_limits(series, spec) # validation before build or output creation
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        input = joinpath(dir,"graph.tsv")
        write_graph_spec(input,series,spec,limits)
        stem = joinpath(dir,"rendered")
        run(Cmd([exe,input,stem]))
        move_plot_outputs(stem, paths)
    end
    return paths
end

function fit_root(data::Union{GraphSeries,Hist1D}, model::Symbol;
                  initial, bounds=nothing, fixed=Int[], range=nothing,
                  objective::Symbol=:chi2, xerrors::Symbol=:effective,
                  root_config=nothing)
    haskey(ROOT_FIT_MODELS,model) || throw(ArgumentError("unsupported fit model: $model"))
    npar = ROOT_FIT_NPAR[model]
    p = Float64.(initial)
    length(p) == npar && all(isfinite,p) || throw(ArgumentError("initial needs $npar finite parameters"))
    lo = fill(-Inf,npar); hi=fill(Inf,npar)
    if bounds !== nothing
        length(bounds) == npar || throw(ArgumentError("one bound pair per parameter is needed"))
        for i in 1:npar
            bounds[i] === nothing && continue
            lo[i],hi[i] = Float64.(bounds[i])
            isfinite(lo[i]) && isfinite(hi[i]) && lo[i] < hi[i] ||
                throw(ArgumentError("bounds must be finite and increasing"))
            lo[i] <= p[i] <= hi[i] || throw(ArgumentError("initial parameter lies outside bounds"))
        end
    end
    all(i -> i isa Integer && 1 <= i <= npar, fixed) && length(unique(fixed)) == length(fixed) ||
        throw(ArgumentError("invalid fixed parameter indices"))
    validate_model_parameters(model,p)
    for i in edge_width_indices(model)
        i in fixed || (isfinite(lo[i]) && lo[i] > 0 && isfinite(hi[i])) ||
            throw(ArgumentError("edge width parameter $i needs finite positive bounds or must be fixed"))
    end
    objective in (:chi2,:poisson,:weighted_poisson) || throw(ArgumentError("invalid objective"))
    xerrors in (:effective,:ignore) || throw(ArgumentError("xerrors must be :effective or :ignore"))
    if data isa GraphSeries
        validate_graph(data)
        objective == :chi2 || throw(ArgumentError("graph fits support :chi2 only"))
        data.kind == :plain || all(>(0),data.eylow) && all(>(0),data.eyhigh) ||
            throw(ArgumentError("graph fit requires strictly positive y errors"))
        default_range = (minimum(data.x),maximum(data.x))
    else
        all(isfinite,data.edges) && all(isfinite,data.sumw) && all(isfinite,data.sumw2) ||
            throw(ArgumentError("histogram data changed to nonfinite values"))
        length(data.edges) == length(data.sumw)+1 == length(data.sumw2)+1 &&
            all(diff(data.edges) .> 0) && all(>=(0),data.sumw2) ||
            throw(ArgumentError("histogram data changed invalidly"))
        if objective != :chi2
            all(>=(0),data.sumw) || throw(ArgumentError("likelihood fit forbids negative bins"))
            if objective == :poisson
                all(x -> x >= 0 && x == round(x),data.sumw) &&
                    all(isapprox.(data.sumw2,data.sumw;atol=1e-10,rtol=1e-10)) ||
                    throw(ArgumentError(":poisson requires unweighted integer counts and sumw2=count"))
            else
                all(i -> data.sumw[i] == 0 || data.sumw2[i] > 0, eachindex(data.sumw)) ||
                    throw(ArgumentError("positive weighted bin requires positive variance"))
            end
        end
        default_range = (first(data.edges),last(data.edges))
    end
    r = range === nothing ? default_range : Tuple(Float64.(range))
    length(r) == 2 && all(isfinite,r) && r[1] < r[2] ||
        throw(ArgumentError("fit range must be two finite increasing values"))
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    mktempdir() do dir
        input,output = joinpath(dir,"fit.tsv"),joinpath(dir,"result.tsv")
        open(input,"w") do io
            println(io,"ROOTPlotsGraph\t1")
            println(io,join(("FIT",data isa GraphSeries ? "graph" : "hist",model,objective,xerrors,
                             r[1],r[2],npar), '\t'))
            for i in 1:npar
                println(io,join(("PARAM",i,p[i],lo[i],hi[i],Int(i in fixed)), '\t'))
            end
            if data isa GraphSeries
                s=data
                println(io,join(("SERIES",1,length(s.x),s.kind,graph_series_draw(s),s.errors,
                                  clean_field(s.label),s.style.color,s.style.marker_style,
                                  s.style.marker_size,s.style.line_width,s.style.line_style), '\t'))
                for i in eachindex(s.x)
                    println(io,join(("POINT",1,i,s.x[i],s.y[i],s.exlow[i],s.exhigh[i],
                                      s.eylow[i],s.eyhigh[i]), '\t'))
                end
            else
                println(io,"HIST\t",length(data.sumw))
                for i in eachindex(data.sumw)
                    println(io,join(("BIN",i,data.edges[i],data.edges[i+1],
                                      data.sumw[i],data.sumw2[i]), '\t'))
                end
            end
            println(io,"END\t1")
        end
        run(Cmd([exe,input,output]))
        return read_fit_result(output,model,r,objective)
    end
end

function read_fit_result(path,model,r,objective)
    lines = readlines(path)
    length(lines) >= 4 && first(lines) == "ROOTPlotsFitResult\t1" &&
        last(lines) == "END" || error("incomplete ROOT fit result")
    meta = split(lines[2],'\t';keepempty=true)
    length(meta) == 10 && meta[1] == "FIT" || error("invalid ROOT fit result")
    result_bool(token) = token == "1" ? true :
                         token == "0" ? false : error("invalid ROOT result boolean")
    npar = ROOT_FIT_NPAR[model]
    length(lines) == 3+2npar || error("invalid fit result record count")
    p=zeros(npar); e=zeros(npar); fixed=falses(npar); cov=zeros(npar,npar)
    for i in 1:npar
        f=split(lines[2+i],'\t')
        length(f)==5 && f[1]=="PARAM" && parse(Int,f[2])==i ||
            error("invalid parameter result")
        p[i]=parse(Float64,f[3]);e[i]=parse(Float64,f[4]);fixed[i]=result_bool(f[5])
    end
    for i in 1:npar
        f=split(lines[2+npar+i],'\t')
        length(f)==npar+2 && f[1]=="COV" && parse(Int,f[2])==i ||
            error("invalid covariance result")
        cov[i,:] = parse.(Float64,f[3:end])
    end
    return FitResult(model,r,p,e,cov,fixed,parse(Int,meta[2]),result_bool(meta[3]),
                     parse(Int,meta[4]),parse(Float64,meta[5]),
                     objective == :chi2 ? :chi2 :
                     objective == :poisson ? :poisson_deviance : :weighted_poisson_objective,
                     parse(Int,meta[6]),parse(Int,meta[7]),parse(Float64,meta[8]),
                     meta[9],result_bool(meta[10]))
end
