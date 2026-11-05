"""
Indexed callback failures from `rootmap`; each captured exception retains its
backtrace. A failed native ROOT process appears as a [`RootBackendError`](@ref)
with that job's native diagnostic (`error.failures[k].second.ex`).
"""
struct RootMapError <: Exception
    failures::Vector{Pair{Int,CapturedException}}
end

function Base.showerror(io::IO, error::RootMapError)
    print(io, "rootmap failed at input indices ", join(first.(error.failures), ", "))
    for (index, failure) in error.failures
        # Multi-line messages (e.g. native stderr) stay indented under their index.
        text = sprint(showerror, failure.ex; context=io)
        print(io, "\n  [", index, "] ", replace(text, "\n" => "\n      "))
    end
end

"""
    rootmap(f, jobs; workers=Threads.nthreads(:default), timeout=nothing)
    rootmap(jobs; workers=4, timeout=120) do job
        # independent fit_root/render_root/histogram calls
    end

Apply `f` using at most `workers` Julia worker tasks on the default thread pool.
Return a vector in collected input order, regardless of completion order. Each
ROOT fit/render still runs in its own native process: this does not enable ROOT
internal multithreading or share a canvas across threads. With one Julia thread,
process/I/O waits can still overlap when `workers > 1`.

Every input is attempted and all workers are joined before returning or throwing
`RootMapError` with failures in input order. Successful output files are retained
on failure; a batch is not transactional. Empty input does not require ROOT.

`timeout` limits each native ROOT backend process started by the callbacks
(seconds; `Inf` for no limit; `nothing` keeps the enclosing `rootmap` batch's
limit, else `ROOTPLOTS_BACKEND_TIMEOUT`, default 600 s). It is fixed when the batch
starts and is passed to the worker tasks through a scoped value, so other
renders and batches are not affected. A job whose backend exceeds it fails like
any other job, with a [`RootBackendError`](@ref) (`timedout=true`, native stderr kept)
in the `RootMapError`; the other jobs continue. Compilation and `root-config` keep
their own limits. The limit covers native processes only: a callback that
blocks in Julia code is not interrupted.

Interrupting the batch (e.g. Ctrl-C while it waits, or an `InterruptException`
in a job) cancels this invocation only: no further job starts, the native
processes started by its jobs (including nested batches) are terminated and
reaped, the workers are joined, and the `InterruptException` is rethrown. Native
processes of concurrent batches or of unrelated calls are not touched. Joining
waits for running callbacks to return; a callback stuck in Julia code keeps the
batch waiting (a second interrupt stops the wait).

Callbacks must be safe to run concurrently. Keep shared data/spec arrays and the
process environment unchanged during execution; use distinct output destinations
for every job (including concurrent calls outside this batch). No deep copies or
duplicate-output preflight are performed. `workers=1` requests sequential jobs.
"""
function rootmap(f, jobs; workers::Integer=Threads.nthreads(:default), timeout=nothing)
    workers >= 1 || throw(ArgumentError("workers must be a positive integer"))
    limit = timeout === nothing ? nothing : timeout_keyword(timeout)
    items = collect(jobs)
    n = length(items)
    n == 0 && return Any[]
    parent = BACKEND_CONTEXT[]
    # Resolved once (invalid environment values fail here, before any job runs).
    limit = something(limit, parent === nothing ? backend_timeout() : parent.timeout)
    context = BackendContext(parent, limit)
    results = Vector{Any}(undef, n)
    failures = Vector{Union{Nothing,CapturedException}}(nothing, n)
    queue = Channel{Int}(n)
    for i in 1:n
        put!(queue, i)
    end
    close(queue)
    worker() = for i in queue
        is_cancelled(context) && break
        try
            results[i] = f(items[i])
        catch exception
            if exception isa InterruptException || is_cancelled(context)
                cancel!(context)          # stop scheduling and this batch's processes
                break
            end
            failures[i] = CapturedException(exception, catch_backtrace())
        end
    end
    tasks = with(BACKEND_CONTEXT => context) do
        [Threads.@spawn :default worker() for _ in 1:min(workers, n)]
    end
    try
        foreach(wait, tasks)
    catch
        cancel!(context)
        for task in tasks
            try wait(task) catch end
        end
        rethrow()
    end
    is_cancelled(context) && throw(InterruptException())
    errors = Pair{Int,CapturedException}[i => failures[i] for i in 1:n if failures[i] !== nothing]
    isempty(errors) || throw(RootMapError(errors))
    return results
end
