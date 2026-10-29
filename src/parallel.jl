"""Indexed callback failures from `rootmap`; each captured exception retains its backtrace."""
struct RootMapError <: Exception
    failures::Vector{Pair{Int,CapturedException}}
end

function Base.showerror(io::IO, error::RootMapError)
    print(io, "rootmap failed at input indices ", join(first.(error.failures), ", "))
    for (index, failure) in error.failures
        print(io, "\n  [", index, "] ")
        showerror(io, failure.ex)
    end
end

"""
    rootmap(f, jobs; workers=Threads.nthreads(:default))
    rootmap(jobs; workers=4) do job
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

Callbacks must be safe to run concurrently. Keep shared data/spec arrays and the
process environment unchanged during execution; use distinct output destinations
for every job (including concurrent calls outside this batch). No deep copies or
duplicate-output preflight are performed. `workers=1` requests sequential jobs.
"""
function rootmap(f, jobs; workers::Integer=Threads.nthreads(:default))
    workers >= 1 || throw(ArgumentError("workers must be a positive integer"))
    items = collect(jobs)
    n = length(items)
    n == 0 && return Any[]
    results = Vector{Any}(undef, n)
    failures = Vector{Union{Nothing,CapturedException}}(nothing, n)
    queue = Channel{Int}(n)
    for i in 1:n
        put!(queue, i)
    end
    close(queue)
    @sync for _ in 1:min(workers, n)
        Threads.@spawn :default begin
            for i in queue
                try
                    results[i] = f(items[i])
                catch exception
                    failures[i] = CapturedException(exception, catch_backtrace())
                end
            end
        end
    end
    errors = Pair{Int,CapturedException}[i => failures[i] for i in 1:n if failures[i] !== nothing]
    isempty(errors) || throw(RootMapError(errors))
    return results
end
