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

"""The one ROOT style named by the comparison panels (all must agree), or nothing."""
function comparison_style(specs)
    styles = [s.style for s in specs if s.style !== nothing]
    length(unique(lowercase.(String.(styles)))) <= 1 ||
        throw(ArgumentError("comparison panels request different ROOT styles; gROOT->SetStyle is process-global"))
    return isempty(styles) ? nothing : first(styles)
end

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
        root_style = comparison_style(specs)
        root_style === nothing || println(io, "ROOTSTYLE\t", String(root_style))
        for (panel_id, (cmp, spec)) in enumerate(zip(cmps, specs))
            n = length(cmp.data.sumw)
            println(io, join(("PANEL", panel_id, root_text(spec.xlabel),
                              root_text(spec.ylabel), root_text(spec.annotation),
                              root_text(spec.lumi), root_text(cmp.uncertainty_label), n), '\t'))
            limits = panel_limits(cmp, spec, ratio_limits)
            println(io, join(("AXIS", panel_id, Int(spec.logx), Int(spec.logy),
                              Int(spec.ratio_logy), limits.ymin, limits.ymax,
                              limits.ratio_min, limits.ratio_max), '\t'))
            write_experiment_label(io, experiment_label(spec; lumi=spec.lumi, annotation=spec.annotation))
            if has_panel_axes(spec)
                validate_panel_axes(spec, n)
                println(io, join(("PCANVAS", panel_id, (scene_value(getfield(spec.canvas, f)) for f in fieldnames(CanvasStyle))...), '\t'))
                for (name, a) in (("x", spec.xaxis), ("y", spec.yaxis), ("ratio", spec.ratio_yaxis))
                    write_axis(io, "PAXIS\t$panel_id", name, a)
                end
            end
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
                println(io, join(("COMP", panel_id, j, root_text(component.label),
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
    fields = split(strip(run_backend(Cmd([exe, "--image-diff", String(reference), String(candidate)]); operation="image comparison", capture=true)), '\t')
    if first(fields) == "SIZE"
        return (; same_size=false, width=parse(Int, fields[2]), height=parse(Int, fields[3]),
                max=Inf, mean=Inf, changed_fraction=1.0)
    end
    first(fields) == "DIFF" && length(fields) == 6 || error("unexpected image comparison output")
    w, h = parse(Int, fields[2]), parse(Int, fields[3])
    return (; same_size=true, width=w, height=h, max=parse(Float64, fields[4]),
            mean=parse(Float64, fields[5]), changed_fraction=parse(Int, fields[6]) / (w * h))
end

"""
    image_colored_fraction(path; root_config=nothing)

Fraction of PNG pixels whose RGB channels differ by more than 32 (internal; used
by tests to check that black-and-white output contains no coloured pixels).
"""
function image_colored_fraction(path::AbstractString; root_config=nothing)
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    fields = split(strip(run_backend(Cmd([exe, "--image-chroma", String(path)]); operation="image colour measurement", capture=true)), '\t')
    first(fields) == "CHROMA" && length(fields) == 4 || error("unexpected image colour output")
    return parse(Int, fields[4]) / (parse(Int, fields[2]) * parse(Int, fields[3]))
end

"""
    image_ink_bounds(path; region=nothing, root_config=nothing)

Bounding box of the dark pixels of a PNG (internal; used by tests to check that
text stays inside a canvas or pad). `region=(x0, y0, x1, y1)` limits the search
to an inclusive pixel rectangle (top-left origin). Returns
`(; width, height, xmin, ymin, xmax, ymax)`, the box fields being `-1` without dark pixels.
"""
function image_ink_bounds(path::AbstractString; region=nothing, root_config=nothing)
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe = graphfit_executable(config)
    args = region === nothing ? String[] : string.(Int.(collect(region)))
    length(args) in (0, 4) || throw(ArgumentError("region must be (x0, y0, x1, y1)"))
    fields = split(strip(run_backend(Cmd([exe, "--image-ink", String(path), args...]); operation="image ink measurement", capture=true)), '\t')
    first(fields) == "INK" && length(fields) == 7 || error("unexpected image ink output")
    v = parse.(Int, fields[2:7])
    return (; width=v[1], height=v[2], xmin=v[3], ymin=v[4], xmax=v[5], ymax=v[6])
end

const ROOT_BUILD_LOCK = ReentrantLock()
# root-config/compiler metadata: separate from the build lock, so cached lookups
# and queries never wait for a compilation.
const ROOT_METADATA_LOCK = ReentrantLock()
const ROOT_METADATA = Dict{Tuple,NamedTuple}()
const COMPILER_VERSIONS = Dict{Tuple,String}()

# Initialized on first build, not at package precompilation time. Scratch owns
# cache lifecycle and keeps all writes outside the (possibly read-only) source.
root_cache_dir() = @get_scratch!("native-backends")

"""Cache key of an executable file: the selected path, its real path, and its device, inode, size and mtime."""
function executable_identity(path::AbstractString)
    real = realpath(path)
    info = stat(real)
    return (String(path), real, info.device, info.inode, info.size, info.mtime)
end

"""
    root_config_metadata(config; force=false)

`root-config --version` and `--cflags --libs` of one installation (internal),
cached per [`executable_identity`](@ref) of the `root-config` file: a replaced or
edited `root-config`, or another path to it, is queried again. The cache does not
notice changes to the ROOT libraries themselves; `force=true` (used by
`build_root(force=true)`) queries again. Queries run under
`ROOTPLOTS_ROOT_CONFIG_TIMEOUT` (default 60 s); one lock prevents duplicate queries.
"""
function root_config_metadata(config::AbstractString; force::Bool=false)
    id = executable_identity(config)
    lock(ROOT_METADATA_LOCK) do
        if !force
            cached = get(ROOT_METADATA, id, nothing)
            cached === nothing || return cached
        end
        real = id[2]
        limit = root_config_timeout()
        version = strip(run_backend(Cmd([real, "--version"]); operation="root-config --version",
                                    capture=true, timeout=limit))
        flags = Base.shell_split(strip(run_backend(Cmd([real, "--cflags", "--libs"]);
                                                   operation="root-config --cflags --libs",
                                                   capture=true, timeout=limit)))
        metadata = (; root_config=real, version=String(version), flags=Vector{String}(flags))
        ROOT_METADATA[id] = metadata
        return metadata
    end
end

"""Forget cached `root-config` and compiler metadata (internal; `build_root(force=true)` refreshes one installation)."""
function clear_backend_cache!()
    lock(ROOT_METADATA_LOCK) do
        empty!(ROOT_METADATA)
        empty!(COMPILER_VERSIONS)
    end
    return nothing
end

"""The C++ compiler selected by `CXX` (default `c++`), resolved on every call."""
function resolve_compiler()
    compiler = Sys.which(get(ENV, "CXX", "c++"))
    compiler === nothing && throw(ArgumentError("C++ compiler not found; set CXX to a compiler executable path (without flags)"))
    return String(compiler)
end

"""First line of `compiler --version`, cached per compiler file identity (provenance only)."""
function compiler_version(compiler::AbstractString)
    id = executable_identity(compiler)
    lock(ROOT_METADATA_LOCK) do
        get!(COMPILER_VERSIONS, id) do
            text = run_backend(Cmd([String(compiler), "--version"]); operation="compiler --version",
                               capture=true, timeout=root_config_timeout())
            String(strip(first(split(text, '\n'; limit=2))))
        end
    end
end

"""
Content/configuration key of one backend build (internal): ROOT version and
flags (cached metadata), compiler path and real path, language standard, and the
SHA-256 of the source and every listed header, hashed again on every call.
"""
function backend_identity(config::AbstractString, source::AbstractString, standard::AbstractString;
                          dependencies=String[])
    metadata = root_config_metadata(config)
    compiler = resolve_compiler()
    source_files = [String(source); String.(dependencies)]
    source_keys = [(realpath(path), bytes2hex(sha256(read(path)))) for path in source_files]
    identity = (metadata.root_config, metadata.version, compiler, realpath(compiler), standard,
                metadata.flags, source_keys)
    return (; key=bytes2hex(sha256(repr(identity))), compiler, metadata)
end

"""Compile privately and atomically publish a content/configuration-keyed backend.

The Julia lock prevents duplicate builds in this process; cache hits do not take
it. Independent Julia processes may compile the same key redundantly, but never
share a staging file or publish a partial binary. No mutable version stamp is
used. The compiler runs under `ROOTPLOTS_COMPILE_TIMEOUT` (default 600 s); a
failure is a [`RootBackendError`](@ref) with the compiler's messages.
"""
function compile_root_backend(config::AbstractString, source::AbstractString,
                              name::AbstractString, standard::AbstractString;
                              force::Bool=false,
                              dependencies=String[],
                              build_dir=root_cache_dir())
    force && root_config_metadata(config; force=true)   # explicit metadata refresh
    id = backend_identity(config, source, standard; dependencies)
    destination = joinpath(build_dir, name, id.key)
    executable = joinpath(destination, name)
    !force && isfile(executable) && return executable
    lock(ROOT_BUILD_LOCK) do
        if force || !isfile(executable)
            limit = compile_timeout()
            mkpath(destination)
            mktempdir(destination; prefix="compile-") do staging
                temporary = joinpath(staging, name)
                run_backend(Cmd(vcat([id.compiler, standard, "-O2", String(source), "-o", temporary], id.metadata.flags));
                            operation="C++ compilation of $name", timeout=limit)
                # Unlike mv(...; force=true), rename replaces without unlinking
                # the live destination first. Existing processes retain their inode.
                Base.Filesystem.rename(temporary, executable)
            end
        end
    end
    return executable
end

"""Compile the local control renderer, reusing a content/configuration-keyed artifact."""
function build_root_backend(config::AbstractString; force::Bool=false)
    source = normpath(joinpath(@__DIR__, "..", "backend", "root", "render.cpp"))
    return compile_root_backend(config, source, "rootplots-render", "-std=c++17"; force,
                                dependencies=[joinpath(dirname(source), header) for header in ("axisstyle.hpp", "experiment.hpp", "outputs.hpp")])
end

"""
    build_root(; root_config=nothing, force=false)

Compile both native backends against a locally installed C++ ROOT using its
`root-config --cflags --libs`. ROOT itself is not built or downloaded.
Resolve ROOT as for `render_root`: explicit `root_config`, `ROOT_CONFIG`,
`ROOTSYS`, then PATH. Return a named tuple containing `root_config`, `version`,
`renderer`, and `graphfit` paths. Compilation/linking failures throw a
[`RootBackendError`](@ref) carrying the compiler's messages.

Unchanged builds are reused; `force=true` queries `root-config` again and recompiles
both backends. Rendering and fitting also build their backend automatically on
first use. The `root-config` output (version and flags) is cached per
`root-config` file (path, real path, inode, size, modification time) for the
session; backend sources and headers are hashed and `CXX` is resolved on every
lookup, so a source edit or another compiler selects another build. Updating
the ROOT libraries behind an unchanged `root-config` is not detected: use
`force=true`. `root-config` runs under `ROOTPLOTS_ROOT_CONFIG_TIMEOUT` (default
60 s), the compiler under `ROOTPLOTS_COMPILE_TIMEOUT` (default 600 s).

Build success does not by itself verify graphics output or scientific
equivalence to a reference plot. Concurrent calls use a locked cache and atomic
binary publication; rendering/fitting processes keep independent ROOT objects
and temporary files.
"""
function build_root(; root_config=nothing, force::Bool=false)
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found; install ROOT or pass root_config to build_root"))
    version = root_config_metadata(config; force).version
    renderer = build_root_backend(config; force)
    graphfit = graphfit_executable(config; force)
    return (; root_config=config, version, renderer, graphfit)
end

"""
    render_root(cmps::AbstractVector{Comparison}, specs::AbstractVector{PanelSpec}, outstem;
                ratio_limits=(0.5, 1.5), style=PlotStyle(), format=:png, provenance=false,
                root_config=nothing)

Render data/model comparisons with ratio panels by the native control renderer
(`rootplots-render`). Julia owns all numeric inputs and no ROOT object enters
the Julia process. `format` (`:png`, `:pdf`, `:both`) selects the files ROOT
writes; `provenance=true` also writes a reproducibility record and returns
`(; png, pdf, request, provenance)` (see the graph method of `render_root`).
"""
function render_root(cmps::AbstractVector{Comparison}, specs::AbstractVector{PanelSpec},
                     outstem::AbstractString; ratio_limits=(0.5, 1.5),
                     style=PlotStyle(), format=:png, provenance::Bool=false, root_config=nothing)
    paths = plot_output_paths(outstem, format)
    length(cmps) == length(specs) && !isempty(cmps) ||
        throw(ArgumentError("one PanelSpec is needed per nonempty comparison"))
    validate_style(style)
    for (cmp, spec) in zip(cmps, specs)
        panel_limits(cmp, spec, ratio_limits)
        validate_panel_axes(spec, length(cmp.data.sumw))
        experiment_label(spec; lumi=spec.lumi, annotation=spec.annotation)
    end
    comparison_style(specs)
    config = find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found; install ROOT or pass root_config to render_root"))
    exe = build_root_backend(config)
    return render_request(paths, outstem, (config, exe); operation="comparison rendering",
                          provenance, name="figure.tsv") do input
        write_root_spec(input, cmps, specs, ratio_limits, style)
    end
end
