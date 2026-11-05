# Central image rendering: one request, one native process, requested formats
# only, and the optional reproducibility record.

"""Relative output stem given to the renderer, which runs in the private temporary directory."""
const RENDERED_STEM = "rendered"

"""
`command` with a relative `TMPDIR` made absolute (against the caller's working
directory) in the child's environment only; unchanged when `TMPDIR` is unset or
absolute. The renderer runs in another directory, where a relative `TMPDIR`
would name a missing directory.
"""
function child_tmpdir(command::Cmd)
    value = get(ENV, "TMPDIR", "")
    isempty(value) || isabspath(value) ? command : addenv(command, "TMPDIR" => abspath(value))
end

"""Native command-line option selecting the requested formats (`--formats=png`, ...)."""
function format_argument(paths)
    formats = String[]
    paths.png === nothing || push!(formats, "png")
    paths.pdf === nothing || push!(formats, "pdf")
    isempty(formats) && throw(ArgumentError("no output format selected"))
    return "--formats=" * join(formats, ",")
end

"""Requested format names of an output-path tuple."""
requested_formats(paths) = [name for (name, path) in (("png", paths.png), ("pdf", paths.pdf)) if path !== nothing]

"""
    render_request(write_request, paths, outstem, backend; operation, provenance=false,
                   name="request.tsv")

Write one request with `write_request(path)` into a private temporary directory,
run the native renderer for the requested formats only, and move the files to
`paths` (internal; used by every image route). `backend` is `(config, exe)` or a
function returning it (called after the request was written, so request
validation errors come first). With `provenance=true` the request and a TOML record
are written next to the images after the rendering succeeded; the return value
then also has `request` and `provenance` paths.

The renderer runs inside the private directory and receives the relative output
stem `"rendered"`, so the file name ROOT stores in the PDF (its `/Title`) is
`rendered.pdf`, never a temporary or home-directory path. The request path stays
absolute; all data are inline in the request. The temporary directory, the request
path and the child's working directory are absolute even when `TMPDIR` is a
relative path; such a `TMPDIR` is also passed to the child as an absolute path
(resolved against the caller's working directory). Julia's own `ENV` is not
changed.
"""
function render_request(write_request, paths, outstem::AbstractString, backend;
                        operation::AbstractString, provenance::Bool=false, name::AbstractString="request.tsv")
    mktempdir() do temporary
        dir = abspath(temporary)            # TMPDIR may be relative
        input = joinpath(dir, name)
        write_request(input)
        config, exe = backend isa Tuple ? backend : backend()
        mkpath(dirname(outstem))
        command = Cmd(Cmd([exe, input, RENDERED_STEM, format_argument(paths)]); dir=dir)
        run_backend(child_tmpdir(command); operation, request=input)
        move_plot_outputs(joinpath(dir, RENDERED_STEM), paths)
        provenance || return paths
        return write_provenance(input, config, exe, outstem, paths; operation)
    end
end

"""Absolute paths in a tool's version line are reduced to their file names."""
strip_version_paths(line::AbstractString) = replace(line, r"(?<![\w.:/])/[^\s()\[\],;]+" => s -> basename(s))

"""Value of the first request row with tag `tag` (field 2), or `""`."""
function request_field(path, tag)
    for line in eachline(path)
        startswith(line, tag * "\t") && return String(split(line, '\t')[2])
    end
    return ""
end

file_sha256(path) = open(io -> bytes2hex(sha256(io)), path)

"""
Write with `writer(io)` to a new temporary file in the destination directory and
return its path. If `writer` (or closing the file) fails, the temporary file is
removed before the exception is rethrown.
"""
function staged_file(writer, destination)
    temporary, io = mktemp(dirname(abspath(destination)); cleanup=false)
    try
        writer(io)
        close(io)
    catch
        close(io)
        rm(temporary; force=true)
        rethrow()
    end
    return temporary
end

"""
Reproducibility record of one successful rendering (internal). The request is
copied byte for byte; the TOML record lists versions, platform, backend key,
style, formats, and SHA-256 digests, with file names but no directory paths,
command lines, host or user names, or environment.

Both files are staged next to their destinations and renamed one after the
other. Each rename is atomic, the pair is not: if the second rename fails, the
request file is already in place (the staged record is removed and the error
rethrown).
"""
function write_provenance(input, config, exe, outstem, paths; operation)
    request_path = outstem * ".request.tsv"
    record_path = outstem * ".provenance.toml"
    request_sha = file_sha256(input)
    metadata = root_config_metadata(config)
    compiler = resolve_compiler()
    outputs = Dict{String,Any}("formats" => requested_formats(paths))
    for (format, path) in (("png", paths.png), ("pdf", paths.pdf))
        path === nothing && continue
        outputs[format] = Dict("file" => basename(path), "sha256" => file_sha256(path))
    end
    protocol = open(readline, input)
    record = Dict{String,Any}(
        "record" => "ROOTPlots provenance 1",
        "rootplots" => Dict("version" => string(pkgversion(@__MODULE__)),
                            "protocol" => replace(protocol, '\t' => ' ')),
        "julia" => Dict("version" => string(VERSION)),
        "root" => Dict("version" => metadata.version),
        "compiler" => Dict("name" => basename(compiler),
                           "version" => strip_version_paths(compiler_version(compiler))),
        "system" => Dict("kernel" => string(Sys.KERNEL), "arch" => string(Sys.ARCH),
                         "machine" => Sys.MACHINE),
        "backend" => Dict("executable" => basename(exe), "key" => basename(dirname(exe)),
                          "operation" => String(operation)),
        "request" => Dict("file" => basename(request_path), "sha256" => request_sha,
                          "kind" => backend_request_kind(input),
                          "style" => request_field(input, "ROOTSTYLE")),
        "outputs" => outputs)
    staged = String[]
    try
        push!(staged, staged_file(io -> write(io, read(input)), request_path))
        push!(staged, staged_file(io -> TOML.print(io, record; sorted=true), record_path))
        Base.Filesystem.rename(staged[1], request_path)
        Base.Filesystem.rename(staged[2], record_path)
    catch
        foreach(path -> rm(path; force=true), staged)
        rethrow()
    end
    return (; png=paths.png, pdf=paths.pdf, request=request_path, provenance=record_path)
end
