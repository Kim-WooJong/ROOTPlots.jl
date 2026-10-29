"""Resolve a local `root-config` executable from an explicit path, ROOT_CONFIG, ROOTSYS, or PATH; return `nothing` when absent without contacting any remote host."""
function find_root_config(explicit::Union{Nothing,AbstractString}=nothing)
    # A deliberately selected installation must not silently fall back to another.
    if explicit !== nothing || haskey(ENV, "ROOT_CONFIG")
        selected = explicit === nothing ? ENV["ROOT_CONFIG"] : String(explicit)
        found = Sys.which(selected)
        return found === nothing ? nothing : String(found)
    end
    candidates = String[]
    haskey(ENV, "ROOTSYS") && push!(candidates, joinpath(ENV["ROOTSYS"], "bin", "root-config"))
    found = Sys.which("root-config")
    found === nothing || push!(candidates, found)
    index = findfirst(isfile, candidates)
    return index === nothing ? nothing : candidates[index]
end

"""Report whether a native local ROOT installation exposes `root-config`; this does not test the compiler, linker, or graphics output."""
root_available(; root_config=nothing) = find_root_config(root_config) !== nothing

"""Reject TSV control characters in caller text; numbers and labels are serialized without shell interpolation or ROOT macro evaluation."""
function clean_field(value::AbstractString)
    any(c -> c in ('\t', '\n', '\r'), value) &&
        throw(ArgumentError("tabs/newlines are not permitted in ROOT backend labels"))
    return value
end

"""Translate the Greek characters used in panel labels into ROOT's TLatex spelling."""
root_latex(value::AbstractString) = replace(String(value), "μ" => "#mu", "ℓ" => "#ell")

"""Write versioned, precomputed plot geometry to a private TSV file for the native renderer. Ratio arithmetic remains Julia-owned; no ROOT files are created."""
function write_root_spec(path::AbstractString, cmps::AbstractVector{Comparison},
                         specs::AbstractVector{PanelSpec}, ratio_limits, style=PlotStyle())
    open(path, "w") do io
        println(io, "ROOTPlots\t2")
        println(io, "SETTINGS\t", Float64(ratio_limits[1]), '\t', Float64(ratio_limits[2]))
        println(io, join(("STYLE", something(style.marker_style, 20),
                          something(style.marker_size, 0.8), something(style.marker_color, "#000000"),
                          something(style.line_width, 1),
                          something(style.stack_fill_style, 1001),
                          something(style.band_fill_style, 3345),
                          something(style.band_color, "#000000"),
                          something(style.legend_columns, 2)), '\t'))
        for (panel_id, (cmp, spec)) in enumerate(zip(cmps, specs))
            n = length(cmp.data.sumw)
            println(io, join(("PANEL", panel_id, clean_field(root_latex(spec.xlabel)),
                              clean_field(root_latex(spec.ylabel)), clean_field(root_latex(spec.annotation)),
                              clean_field(root_latex(spec.lumi)), clean_field(root_latex(cmp.uncertainty_label)), n), '\t'))
            limits = panel_limits(cmp, spec, ratio_limits)
            println(io, join(("AXIS", panel_id, Int(spec.logx), Int(spec.logy),
                              Int(spec.ratio_logy), limits.ymin, limits.ymax,
                              limits.ratio_min, limits.ratio_max), '\t'))
            for i in 1:n
                valid = cmp.valid[i]
                println(io, join(("BIN", panel_id, i, cmp.data.edges[i], cmp.data.edges[i+1],
                                  cmp.data.sumw[i], cmp.data_errors.low[i],
                                  cmp.data_errors.high[i], cmp.prediction[i],
                                  cmp.prediction_errors.low[i], cmp.prediction_errors.high[i],
                                  valid ? cmp.ratio[i] : 0.0,
                                  valid ? cmp.ratio_data_errors.low[i] : 0.0,
                                  valid ? cmp.ratio_data_errors.high[i] : 0.0,
                                  valid ? cmp.ratio_prediction_errors.low[i] : 0.0,
                                  valid ? cmp.ratio_prediction_errors.high[i] : 0.0,
                                  Int(valid)), '\t'))
            end
            for (j, component) in enumerate(cmp.components)
                println(io, join(("COMP", panel_id, j, clean_field(root_latex(component.label)),
                                  component.color, component.hist.sumw...), '\t'))
            end
        end
    end
    return path
end

"""
    image_difference(reference, candidate; root_config=nothing)

Compare two PNG files with ROOT's image reader (internal; used by the visual
regression tests). Returns `(; same_size, width, height, max, mean, changed_fraction)`
where `max`/`mean` are RGB channel differences on a 0–255 scale and a pixel is
counted as changed when any channel differs by more than 32.
"""
function image_difference(reference::AbstractString, candidate::AbstractString; root_config=nothing)
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    fields = split(strip(read(Cmd([exe, "--image-diff", String(reference), String(candidate)]), String)), '\t')
    if first(fields) == "SIZE"
        return (; same_size=false, width=parse(Int, fields[2]), height=parse(Int, fields[3]),
                max=Inf, mean=Inf, changed_fraction=1.0)
    end
    first(fields) == "DIFF" && length(fields) == 6 || error("unexpected image comparison output")
    w, h = parse(Int, fields[2]), parse(Int, fields[3])
    return (; same_size=true, width=w, height=h, max=parse(Float64, fields[4]),
            mean=parse(Float64, fields[5]), changed_fraction=parse(Int, fields[6]) / (w * h))
end

const ROOT_BUILD_LOCK = ReentrantLock()

# Initialized on first build, not at package precompilation time. Scratch owns
# cache lifecycle and keeps all writes outside the (possibly read-only) source.
root_cache_dir() = @get_scratch!("native-backends")

"""Compile privately and atomically publish a content/configuration-keyed backend.

The Julia lock prevents duplicate builds in this process. Independent Julia
processes may compile the same key redundantly, but never share a staging file
or publish a partial binary. No mutable version stamp is used.
"""
function compile_root_backend(config::AbstractString, source::AbstractString,
                              name::AbstractString, standard::AbstractString;
                              force::Bool=false,
                              dependencies=String[],
                              build_dir=root_cache_dir())
    lock(ROOT_BUILD_LOCK) do
        config = realpath(config)
        compiler = Sys.which(get(ENV, "CXX", "c++"))
        compiler === nothing && throw(ArgumentError("C++ compiler not found; set CXX to a compiler executable path (without flags)"))
        version = strip(read(Cmd([config, "--version"]), String))
        flags = Base.shell_split(strip(read(Cmd([config, "--cflags", "--libs"]), String)))
        source_files = [String(source); String.(dependencies)]
        source_keys = [(realpath(path), bytes2hex(sha256(read(path)))) for path in source_files]
        identity = (config, version, compiler, realpath(compiler), standard, flags, source_keys)
        key = bytes2hex(sha256(repr(identity)))
        destination = joinpath(build_dir, name, key)
        executable = joinpath(destination, name)
        if force || !isfile(executable)
            mkpath(destination)
            mktempdir(destination; prefix="compile-") do staging
                temporary = joinpath(staging, name)
                run(Cmd(vcat([compiler, standard, "-O2", String(source), "-o", temporary], flags)))
                # Unlike mv(...; force=true), rename replaces without unlinking
                # the live destination first. Existing processes retain their inode.
                Base.Filesystem.rename(temporary, executable)
            end
        end
        return executable
    end
end

"""Compile the local control renderer, reusing a content/configuration-keyed artifact."""
function build_root_backend(config::AbstractString; force::Bool=false)
    source = normpath(joinpath(@__DIR__, "..", "backend", "root", "render.cpp"))
    return compile_root_backend(config, source, "rootplots-render", "-std=c++17"; force)
end

"""
    build_root(; root_config=nothing, force=false)

Compile both native backends against a locally installed C++ ROOT using its
`root-config --cflags --libs`. ROOT itself is not built or downloaded.
Resolve ROOT as for `render_root`: explicit `root_config`, `ROOT_CONFIG`,
`ROOTSYS`, then PATH. Return a named tuple containing `root_config`, `version`,
`renderer`, and `graphfit` paths. Compilation/linking failures propagate.

Unchanged builds are reused; `force=true` recompiles both backends. Rendering
and fitting also build their backend automatically on first use. Build success
does not by itself verify graphics output or scientific equivalence to a
reference plot. Concurrent calls use a locked cache and atomic binary publication;
rendering/fitting processes keep independent ROOT objects and temporary files.
"""
function build_root(; root_config=nothing, force::Bool=false)
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found; install ROOT or pass root_config to build_root"))
    version = strip(read(Cmd([config, "--version"]), String))
    renderer = build_root_backend(config; force)
    graphfit = graphfit_executable(config; force)
    return (; root_config=config, version, renderer, graphfit)
end

"""Render selected PNG/PDF output with installed native ROOT through a compiled, batch-mode process. Julia owns all numeric inputs and no ROOT object enters the Julia process."""
function render_root(cmps::AbstractVector{Comparison}, specs::AbstractVector{PanelSpec},
                     outstem::AbstractString; ratio_limits=(0.5, 1.5),
                     style=PlotStyle(), format=:png, root_config=nothing)
    paths = plot_output_paths(outstem, format)
    length(cmps) == length(specs) && !isempty(cmps) ||
        throw(ArgumentError("one PanelSpec is needed per nonempty comparison"))
    validate_style(style)
    for (cmp, spec) in zip(cmps, specs)
        panel_limits(cmp, spec, ratio_limits)
    end
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found; install ROOT or pass root_config to render_root"))
    exe = build_root_backend(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        specfile = write_root_spec(joinpath(dir, "figure.tsv"), cmps, specs, ratio_limits, style)
        temporary_stem = joinpath(dir, "rendered")
        run(Cmd([exe, specfile, temporary_stem]))
        move_plot_outputs(temporary_stem, paths)
    end
    return paths
end
