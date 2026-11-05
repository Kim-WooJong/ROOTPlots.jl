"""
    RootBackendError <: Exception

A native process run by the package failed: a ROOT backend
(`rootplots-graphfit` or `rootplots-render`), the C++ compiler that builds them, or
`root-config`. Thrown by every operation that runs such a process (`render_root`,
`fit_root`, `profile_stats`, the special-graph point helpers, `root_styles`, the
image helpers, `build_root` and the automatic builds) instead of a bare
`ProcessFailedException`, so the native diagnostic survives the deletion of the
temporary request file.

Fields:

- `operation`: what the package was doing (e.g. `"graph rendering"`).
- `request`: the first record of the request file (e.g. `"DRAWFIT"`), `""`
  when the command has no request file.
- `executable`: the executable's file name.
- `exitcode`: the process exit code.
- `termsignal`: the terminating signal number, `0` for a normal exit.
- `stderr`: the native standard error text (ROOT and backend messages). Very
  long output keeps its first and last 32 KiB around an omission marker.
- `timedout`: `true` when the package stopped the process because it exceeded its
  time limit (see "Time limits" below); `termsignal` then
  names the signal that ended it (normally 15, or 9 when it ignored SIGTERM).
- `timeout`: the limit in seconds that applied (`Inf` without a limit).

Standard error is captured while the process runs (large output cannot block
it). After a successful run it is written to Julia's `stderr`, so ROOT warnings
and `ROOTPLOTS_ROOT_VERBOSE=1` messages stay visible.

# Time limits

Every native process runs under a wall-clock limit, read when it starts:

| process | environment variable | default |
|:--|:--|:--|
| ROOT backends (render, fit, statistics, images, styles) | `ROOTPLOTS_BACKEND_TIMEOUT` | 600 s |
| `root-config` queries | `ROOTPLOTS_ROOT_CONFIG_TIMEOUT` | 60 s |
| C++ compilation of a backend | `ROOTPLOTS_COMPILE_TIMEOUT` | 600 s |

A value is a positive number of seconds; `0` or `none` disables the limit; an
empty or unset variable means the default. Anything else is an `ArgumentError`
before a process is started. `rootmap(...; timeout=t)` replaces the backend
limit for the native backend processes of that batch.

On POSIX systems (macOS, Linux) each process is started as the leader of its
own process group (session). When the limit is exceeded, or the waiting Julia
task is interrupted, the package sends SIGTERM to that group, waits up to 2 s,
sends SIGKILL to the group if anything is left, and reaps the process; it never
signals the Julia process's own group. On Windows only the direct child is
terminated (descendants it started may survive). The limit applies to native
processes only: Julia code (including `rootmap` callbacks) cannot be timed out.
"""
struct RootBackendError <: Exception
    operation::String
    request::String
    executable::String
    exitcode::Int
    termsignal::Int
    stderr::String
    timedout::Bool
    timeout::Float64
end

# The six-field constructor of earlier versions: a failure without a time limit event.
RootBackendError(operation, request, executable, exitcode, termsignal, stderr) =
    RootBackendError(operation, request, executable, exitcode, termsignal, stderr, false, Inf)

const BACKEND_STDERR_KEEP = 32 * 1024
# Successful runs forward up to this much standard error verbatim (head and tail beyond it).
const BACKEND_STDERR_FORWARD = 64 * 1024 * 1024
const SIGNAL_NAMES = Dict(1 => "SIGHUP", 2 => "SIGINT", 3 => "SIGQUIT", 4 => "SIGILL",
                          6 => "SIGABRT", 8 => "SIGFPE", 9 => "SIGKILL", 11 => "SIGSEGV",
                          13 => "SIGPIPE", 14 => "SIGALRM", 15 => "SIGTERM")

function Base.showerror(io::IO, e::RootBackendError)
    print(io, "RootBackendError: native ROOT backend failed during ", e.operation)
    isempty(e.request) || print(io, " (request ", e.request, ")")
    if e.timedout
        print(io, ": timed out after ", format_seconds(e.timeout), " s")
        if e.termsignal != 0
            name = get(SIGNAL_NAMES, e.termsignal, "")
            print(io, "; stopped by signal ", e.termsignal, isempty(name) ? "" : " ($name)")
        end
    elseif e.termsignal != 0
        name = get(SIGNAL_NAMES, e.termsignal, "")
        print(io, ": terminated by signal ", e.termsignal, isempty(name) ? "" : " ($name)")
    else
        print(io, ": exit code ", e.exitcode)
    end
    print(io, "\n  executable: ", e.executable)
    text = rstrip(e.stderr)
    if isempty(text)
        print(io, "\n  native stderr: (empty)")
    else
        print(io, "\n  native stderr:")
        for line in split(text, '\n')
            print(io, "\n    ", line)
        end
    end
end

format_seconds(t) = isinteger(t) && abs(t) < 1e15 ? string(Int(t)) : string(t)

"""Keep the head and tail of very long native output."""
function bounded_backend_text(bytes::Vector{UInt8})
    n = length(bytes)
    n <= 2BACKEND_STDERR_KEEP && return String(bytes)
    head = String(bytes[1:BACKEND_STDERR_KEEP])
    tail = String(bytes[n-BACKEND_STDERR_KEEP+1:n])
    return string(head, "\n[... ", n - 2BACKEND_STDERR_KEEP, " bytes of native output omitted ...]\n", tail)
end

"""First record tag of a protocol request file (`""` when unreadable)."""
function backend_request_kind(path)
    path === nothing && return ""
    try
        open(path) do io
            readline(io)                       # format header
            return String(first(split(readline(io), '\t')))
        end
    catch
        return ""
    end
end

"""Reject a request that would hand non-ASCII text to ROOT (TLatex would draw mojibake)."""
function check_ascii_request(path)
    path === nothing && return nothing
    for (number, line) in enumerate(eachline(path))
        isascii(line) && continue
        bad = first(c for c in line if !isascii(c))
        throw(ArgumentError("internal error: untranslated non-ASCII text '$bad' reached the ROOT request " *
                            "(line $number); every label must pass through ROOTPlots.root_latex"))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Time limits

const DEFAULT_BACKEND_TIMEOUT = 600.0
const DEFAULT_ROOT_CONFIG_TIMEOUT = 60.0
const DEFAULT_COMPILE_TIMEOUT = 600.0
# SIGTERM -> SIGKILL grace, and the wait for output streams after the child exits.
const BACKEND_TERM_GRACE = 2.0
const BACKEND_DRAIN_GRACE = 2.0

"""Parse a time-limit environment variable: seconds > 0, or `Inf` for `0`/`none`."""
function timeout_setting(name::AbstractString, default::Float64)
    raw = get(ENV, name, "")
    text = lowercase(strip(raw))
    isempty(text) && return default
    text == "none" && return Inf
    value = tryparse(Float64, text)
    (value === nothing || !isfinite(value) || value < 0) &&
        throw(ArgumentError("$name must be a positive number of seconds, 0, or none (got $(repr(raw)))"))
    return value == 0 ? Inf : value
end

"""Validate a `timeout` keyword: positive seconds, `Inf` for no limit."""
function timeout_keyword(value)
    value isa Real && !isnan(value) && value > 0 ||
        throw(ArgumentError("timeout must be a positive number of seconds (Inf for no limit), got $(repr(value))"))
    return Float64(value)
end

root_config_timeout() = timeout_setting("ROOTPLOTS_ROOT_CONFIG_TIMEOUT", DEFAULT_ROOT_CONFIG_TIMEOUT)
compile_timeout() = timeout_setting("ROOTPLOTS_COMPILE_TIMEOUT", DEFAULT_COMPILE_TIMEOUT)

"""Limit for ROOT backend processes: the enclosing `rootmap` batch's, else the environment's."""
function backend_timeout()
    context = BACKEND_CONTEXT[]
    context === nothing || return context.timeout
    return timeout_setting("ROOTPLOTS_BACKEND_TIMEOUT", DEFAULT_BACKEND_TIMEOUT)
end

# ---------------------------------------------------------------------------
# Process ownership

"""
One native child process started by `run_backend` (internal). `group` is true when
the child was verified to lead its own process group (POSIX), which is then
signalled as a whole. The first stop request starts the single termination task.
"""
mutable struct NativeProcess
    process::Base.Process
    pid::Int32
    group::Bool
    lock::ReentrantLock
    reason::Symbol                  # :none, :timeout, :cancelled, :interrupted, :failed, :stray
    stopper::Union{Nothing,Task}
end

"""
Cancellation scope of one `rootmap` invocation (internal), visible to its worker
tasks through the scoped value `BACKEND_CONTEXT`. `processes` are the live native
processes started inside the scope (including nested batches); `cancel!` stops
only those, and no new process starts in a cancelled scope or its descendants.
"""
mutable struct BackendContext
    lock::ReentrantLock
    parent::Union{Nothing,BackendContext}
    timeout::Float64
    cancelled::Bool
    processes::Set{NativeProcess}
end

BackendContext(parent, timeout) = BackendContext(ReentrantLock(), parent, Float64(timeout), false, Set{NativeProcess}())

const BACKEND_CONTEXT = ScopedValue{Union{Nothing,BackendContext}}(nothing)

"""Contexts from the outermost to `context` (locks are always taken in this order)."""
function context_chain(context)
    chain = BackendContext[]
    while context !== nothing
        pushfirst!(chain, context)
        context = context.parent
    end
    return chain
end

is_cancelled(context) = any(c -> c.cancelled, context_chain(context))

"""Cancel a scope: no new process starts in it, and its live processes are stopped."""
function cancel!(context::BackendContext)
    owned = lock(context.lock) do
        context.cancelled = true
        collect(context.processes)
    end
    foreach(h -> request_stop!(h, :cancelled), owned)
    return nothing
end

const SIGTERM_NUMBER = Cint(15)
const SIGKILL_NUMBER = Cint(9)

posix_group_supported() = Sys.isunix()

"""True when `pid` leads its own process group, which is not this process's group."""
function owns_process_group(pid)
    posix_group_supported() || return false
    group = ccall(:getpgid, Cint, (Cint,), pid)
    return group == pid && group != ccall(:getpgrp, Cint, ())
end

"""Does any process of group `pid` remain (signal 0 probe)?"""
group_alive(pid) = ccall(:killpg, Cint, (Cint, Cint), pid, 0) == 0

"""Send `signal` to the owned process group, or to the direct child only."""
function signal_native(h::NativeProcess, signal::Cint)
    if h.group
        # The group id stays reserved while any member exists; after the leader
        # was reaped, signal only when members remain.
        process_exited(h.process) && !group_alive(h.pid) && return
        ccall(:killpg, Cint, (Cint, Cint), h.pid, signal)
    else
        process_exited(h.process) || kill(h.process, signal)
    end
    return nothing
end

"""SIGTERM, bounded grace, then SIGKILL for whatever is left of the owned process (group)."""
function terminate_native(h::NativeProcess)
    signal_native(h, SIGTERM_NUMBER)
    deadline = time() + BACKEND_TERM_GRACE
    while time() < deadline
        process_exited(h.process) && !(h.group && group_alive(h.pid)) && return nothing
        sleep(0.02)
    end
    signal_native(h, SIGKILL_NUMBER)
    # Bounded wait until the killed group has gone (its members are reaped by
    # Julia or, after reparenting, by the system).
    deadline = time() + BACKEND_TERM_GRACE
    while time() < deadline
        process_exited(h.process) && !(h.group && group_alive(h.pid)) && return nothing
        sleep(0.02)
    end
    return nothing
end

"""Start the termination of `h` once; later requests keep the first reason."""
function request_stop!(h::NativeProcess, reason::Symbol)
    lock(h.lock) do
        h.stopper === nothing || return
        h.reason = reason
        h.stopper = errormonitor(@async terminate_native(h))
    end
    return nothing
end

wait_stopper(h::NativeProcess) = (task = lock(() -> h.stopper, h.lock); task === nothing || wait(task); nothing)

"""
Stop `h` after its waiting task failed (interrupt, cancellation, or an error) and
wait until the stopper finished and the child was reaped. Further interrupts
during this bounded wait are ignored, so the process is never released from its
scopes while it (or its group) may still run.
"""
function stop_after_failure(h::NativeProcess, exception)
    request_stop!(h, exception isa InterruptException ? :interrupted : :failed)
    while true
        try
            wait_stopper(h)
            wait(h.process)
            return nothing
        catch e
            e isa InterruptException || rethrow()
        end
    end
end

"""Start `command` and register it in every enclosing scope (refused once cancelled)."""
function launch_native(command, context)
    chain = context_chain(context)
    foreach(c -> lock(c.lock), chain)
    try
        any(c -> c.cancelled, chain) && throw(InterruptException())
        process = run(command; wait=false)
        pid = getpid(process)
        h = NativeProcess(process, pid, owns_process_group(pid), ReentrantLock(), :none, nothing)
        foreach(c -> push!(c.processes, h), chain)
        return h
    finally
        foreach(c -> unlock(c.lock), reverse(chain))
    end
end

function release_native(h::NativeProcess, context)
    for c in context_chain(context)
        lock(() -> delete!(c.processes, h), c.lock)
    end
end

# ---------------------------------------------------------------------------
# Output capture

"""Standard error kept verbatim up to `BACKEND_STDERR_FORWARD`, then only its tail."""
mutable struct NativeText
    head::Vector{UInt8}
    tail::Vector{UInt8}
    total::Int
end
NativeText() = NativeText(UInt8[], UInt8[], 0)

function append_text!(t::NativeText, chunk::Vector{UInt8})
    t.total += length(chunk)
    room = BACKEND_STDERR_FORWARD - length(t.head)
    if room > 0
        n = min(room, length(chunk))
        append!(t.head, view(chunk, 1:n))
        chunk = chunk[n+1:end]
    end
    isempty(chunk) && return t
    append!(t.tail, chunk)
    excess = length(t.tail) - BACKEND_STDERR_KEEP
    excess > 0 && deleteat!(t.tail, 1:excess)
    return t
end

"""All text when complete; otherwise the head, an omission marker, and the tail."""
function full_text(t::NativeText)
    t.total == length(t.head) && return t.head
    omitted = t.total - length(t.head) - length(t.tail)
    return vcat(t.head, Vector{UInt8}("\n[... $omitted bytes of native output omitted ...]\n"), t.tail)
end

"""Error text: head and tail of at most `BACKEND_STDERR_KEEP` bytes each."""
function error_text(t::NativeText)
    t.total == length(t.head) && return bounded_backend_text(t.head)
    head = String(t.head[1:min(end, BACKEND_STDERR_KEEP)])
    omitted = t.total - length(head) - length(t.tail)
    return string(head, "\n[... ", omitted, " bytes of native output omitted ...]\n", String(copy(t.tail)))
end

function read_native_text(pipe)
    text = NativeText()
    try
        while !eof(pipe)
            append_text!(text, readavailable(pipe))
        end
    catch exception
        exception isa Base.IOError || exception isa EOFError || rethrow()
    end
    return text
end

function read_native_bytes(pipe)
    buffer = IOBuffer()
    try
        write(buffer, pipe)
    catch exception
        exception isa Base.IOError || exception isa EOFError || rethrow()
    end
    return take!(buffer)
end

# ---------------------------------------------------------------------------

"""
    run_backend(cmd; operation, request=nothing, capture=false, timeout=nothing)

Run one native process (internal; the only place that starts ROOT backends,
the compiler, and `root-config`). Standard error is captured concurrently; on
failure a [`RootBackendError`](@ref) with the native text, the exit code or
signal, and `operation` is thrown; on success the text is forwarded to Julia's
`stderr`. `capture=true` returns standard output as a `String` (otherwise it is
inherited and `nothing` is returned); standard input is `devnull`. `request`
names the protocol file: it gives the error context and is checked to contain
only ASCII.

`timeout` (seconds, `Inf` = none; default `backend_timeout()`) is enforced
by stopping the process (group): `RootBackendError` with `timedout=true`. An
interrupt of the waiting task, or the cancellation of the enclosing `rootmap`
scope, stops it the same way and throws `InterruptException`. Descendants that
still hold the output streams 2 s after the child exited are stopped as well.
"""
function run_backend(cmd::Cmd; operation::AbstractString, request=nothing, capture::Bool=false,
                     timeout::Union{Nothing,Real}=nothing)
    limit = timeout === nothing ? backend_timeout() : Float64(timeout)
    limit > 0 || throw(ArgumentError("timeout must be positive"))
    check_ascii_request(request)
    context = BACKEND_CONTEXT[]
    errpipe = Pipe()
    outpipe = capture ? Pipe() : nothing
    child = Cmd(cmd; ignorestatus=true, detach=true)
    command = capture ? pipeline(child; stdin=devnull, stdout=outpipe, stderr=errpipe) :
                        pipeline(child; stdin=devnull, stderr=errpipe)
    h = launch_native(command, context)
    timer = nothing
    drain = nothing
    try
        close(errpipe.in)
        capture && close(outpipe.in)
        errtask = @async read_native_text(errpipe)
        outtask = capture ? (@async read_native_bytes(outpipe)) : nothing
        isfinite(limit) && (timer = Timer(_ -> request_stop!(h, :timeout), limit))
        # The whole lifecycle (waiting for the child, then draining its streams)
        # is covered: an interrupt or error at any point stops the owned process
        # (group) and waits for it before the process leaves its scopes. The
        # original exception is rethrown unchanged.
        text, output = try
            wait(h.process)
            timer === nothing || close(timer)
            # Descendants may keep the streams open after the child exited: stop the
            # owned group, then give up on the streams.
            drain = Timer(BACKEND_DRAIN_GRACE) do _
                errormonitor(@async begin
                    request_stop!(h, :stray)
                    wait_stopper(h)
                    sleep(BACKEND_DRAIN_GRACE)
                    close(errpipe)
                    outpipe === nothing || close(outpipe)
                end)
            end
            drained = (fetch(errtask), outtask === nothing ? nothing : fetch(outtask))
            close(drain)
            wait_stopper(h)
            drained
        catch exception
            stop_after_failure(h, exception)
            rethrow()
        end
        process = h.process
        if h.reason == :cancelled
            throw(InterruptException())
        end
        failed = process.exitcode != 0 || process.termsignal != 0
        timedout = h.reason == :timeout && failed
        if failed
            throw(RootBackendError(String(operation), backend_request_kind(request),
                                   basename(first(cmd.exec)), Int(process.exitcode),
                                   Int(process.termsignal), error_text(text), timedout, limit))
        end
        h.reason == :stray &&
            @warn "native process $(basename(first(cmd.exec))) left descendants holding its output streams; they were terminated" operation
        bytes = full_text(text)
        isempty(bytes) || write(stderr, bytes)
        return capture ? String(output) : nothing
    finally
        timer === nothing || close(timer)
        drain === nothing || close(drain)
        close(errpipe)
        outpipe === nothing || close(outpipe)
        release_native(h, context)
    end
end
