"""Sample a Julia callable as a display-only ROOT TGraph, not an analytic TF1.
`n` fixes the number of equally spaced x samples; no fit or adaptive integration occurs.
"""
function function_curve(f, xrange; n::Int=300, label="", style=GraphStyle(), draw=:line)
    2 <= n <= 1_000_000 || throw(ArgumentError("function curve needs 2:1,000,000 samples"))
    r = Tuple(Float64.(xrange))
    length(r)==2 && all(isfinite,r) && r[1]<r[2] ||
        throw(ArgumentError("function curve range must be finite and increasing"))
    x = collect(range(r...;length=n))
    raw_y = f.(x)
    all(v -> v isa Real && isfinite(v),raw_y) ||
        throw(ArgumentError("function curve must return finite real values"))
    return GraphSeries(x,Float64.(raw_y);label,style,draw,errors=:none)
end

"""A genuine ROOT TF1 from the package's predefined models (not arbitrary C++).
`npx` controls ROOT's drawing resolution; parameters follow `fit_root` ordering.
"""
struct ModelCurve
    model::Symbol
    parameters::Vector{Float64}
    range::Tuple{Float64,Float64}
    npx::Int
    label::String
    style::GraphStyle
end

function ModelCurve(model::Symbol, parameters::AbstractVector{<:Real};
                    range, npx::Int=300, label="", style=GraphStyle(color="#cc3311"))
    r=Tuple(Float64.(range))
    length(r)==2 && all(isfinite,r) && r[1]<r[2] ||
        throw(ArgumentError("model curve range must be finite and increasing"))
    curve=ModelCurve(model,Float64.(parameters),r,npx,String(label),style)
    validate_model_curve(curve)
    return curve
end

function validate_model_curve(curve::ModelCurve)
    validate_model_parameters(curve.model,curve.parameters)
    curve.model==:gaus && curve.parameters[3]==0 &&
        throw(ArgumentError("Gaussian model curve width cannot be zero"))
    all(isfinite,curve.range) && curve.range[1]<curve.range[2] && 2<=curve.npx<=1_000_000 ||
        throw(ArgumentError("invalid model curve range/resolution"))
    scene_text(curve.label)
    check_graph_style(curve.style)
    return curve
end

function model_curve_limits(curve::ModelCurve,spec::GraphSpec)
    validate_model_curve(curve)
    # Only frame bounds use Julia samples; the displayed object below is TF1.
    samples=function_curve(x -> model_value(curve.model,curve.parameters,x),curve.range;
                           n=curve.npx,label=curve.label,style=curve.style)
    return graph_limits([samples],spec)
end

function write_function_spec(path,curve::ModelCurve,spec::GraphSpec,limits)
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        println(io,join(("DRAWFUNCTION",spec.size...,Int(spec.logx),Int(spec.logy),
                         limits[1]...,limits[2]...,String(spec.legend)),'\t'))
        println(io,join(("LABELS",scene_text(spec.title),scene_text(spec.xlabel),
                         scene_text(spec.ylabel)),'\t'))
        write_scene(io,spec)
        println(io,join(("FUNCTION",curve.model,curve.range...,curve.npx,
                         length(curve.parameters),scene_text(curve.label),
                         curve.style.color,curve.style.line_width,curve.style.line_style,
                         curve.parameters...),'\t'))
        println(io,"END\t1")
    end
    return path
end
write_function_spec(path,curve::ModelCurve,spec::GraphSpec) =
    write_function_spec(path,curve,spec,model_curve_limits(curve,spec))

function render_root(curve::ModelCurve,spec::GraphSpec,outstem::AbstractString;
                     format=:png,root_config=nothing)
    paths=plot_output_paths(outstem,format)
    limits=model_curve_limits(curve,spec)
    config=find_root_config(root_config)
    config===nothing && throw(ArgumentError("local native ROOT not found"))
    exe=graphfit_executable(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        input=write_function_spec(joinpath(dir,"function.tsv"),curve,spec,limits)
        stem=joinpath(dir,"rendered")
        run(Cmd([exe,input,stem]))
        move_plot_outputs(stem,paths)
    end
    return paths
end
