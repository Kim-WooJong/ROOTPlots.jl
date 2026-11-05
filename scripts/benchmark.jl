# Warmed micro-benchmarks of ROOTPlots' own overheads (no BenchmarkTools).
#
#   julia --threads=4 --project=. scripts/benchmark.jl [output.tsv]
#
# Environment (all optional):
#   ROOTPLOTS_BENCH_REPS    repetitions per case after one warm-up run (default 7)
#   ROOTPLOTS_BENCH_EVENTS  events for the histogram-fill case (default 1_000_000;
#                           e.g. 10_000_000 for a larger run)
#   ROOTPLOTS_BENCH_CASES   comma-separated subset of: lookup, small_hist, large_heat,
#                           fill, rootmap (default all)
#
# Every case is run once to warm up (compilation, first ROOT start), then timed
# ROOTPLOTS_BENCH_REPS times; the median and interquartile range (IQR) are
# reported in milliseconds. Only APIs present since 0.2 are used, so the same
# script can time an older source tree (run it with --project=<that tree>).
# The report names Julia/ROOT versions, OS kernel, architecture and threads,
# never local paths. Process start-up and file-system noise dominate small
# cases: compare runs only from the same machine and session conditions.
using ROOTPlots
using Statistics: median, quantile
using Printf: @sprintf

const REPS = parse(Int, get(ENV, "ROOTPLOTS_BENCH_REPS", "7"))
const EVENTS = parse(Int, replace(get(ENV, "ROOTPLOTS_BENCH_EVENTS", "1000000"), "_" => ""))
const CASES = split(get(ENV, "ROOTPLOTS_BENCH_CASES", "lookup,small_hist,large_heat,fill,rootmap"), ',')
REPS >= 3 || error("ROOTPLOTS_BENCH_REPS must be at least 3")

# Warm up once, then time f REPS times (seconds).
function measure(f; reps=REPS)
    f()
    return [(@elapsed f()) for _ in 1:reps]
end

rows = Tuple{String,String,Float64,Float64,Int}[]
function report!(case, variant, seconds; per=1)
    ms = 1000 .* seconds ./ per
    q1, q3 = quantile(ms, 0.25), quantile(ms, 0.75)
    push!(rows, (case, variant, median(ms), q3 - q1, length(ms)))
    println(rpad(case, 12), rpad(variant, 30), @sprintf("median %10.3f ms   IQR %9.3f ms   n=%d", median(ms), q3 - q1, length(ms)))
end

built = build_root()
println("ROOTPlots benchmark")
println("  ROOTPlots ", pkgversion(ROOTPlots), ", Julia ", VERSION, ", ROOT ", built.version)
println("  OS ", Sys.KERNEL, ", arch ", Sys.ARCH, ", Julia threads ", Threads.nthreads(:default),
        ", CPU threads ", Sys.CPU_THREADS)
println("  repetitions ", REPS, " (after one warm-up), histogram events ", EVENTS)

mktempdir() do dir
    if "lookup" in CASES
        # Cache-hit lookup of a compiled backend (what every render/fit does first).
        config = built.root_config
        n = 50
        report!("lookup", "graphfit (per call)", measure(() -> for _ in 1:n
            ROOTPlots.graphfit_executable(config) end); per=n)
        report!("lookup", "render (per call)", measure(() -> for _ in 1:n
            ROOTPlots.build_root_backend(config) end); per=n)
    end
    if "small_hist" in CASES
        h = Hist1D(collect(0.0:1.0:50.0), collect(1.0:50.0), collect(1.0:50.0))
        for format in (:png, :pdf, :both)
            report!("small_hist", "TH1D E1 $format", measure(() ->
                render_root(h, GraphSpec(xlabel="x", ylabel="Events"), joinpath(dir, "small-$format"); draw="E1", format)))
        end
    end
    if "large_heat" in CASES
        nx = ny = 400
        values = [sin(i / 25) * cos(j / 31) + 2 for i in 1:nx, j in 1:ny]
        h2 = Hist2D(collect(0.0:nx), collect(0.0:ny), values, fill(0.01, nx, ny))
        for format in (:png, :pdf, :both)
            report!("large_heat", "TH2D 400x400 COLZ $format", measure(() ->
                render_root(h2, Graph2DSpec(), joinpath(dir, "heat-$format"); draw="COLZ", format)))
        end
    end
    if "fill" in CASES
        x = [mod(0.6180339887 * i, 1.0) * 100 for i in 1:EVENTS]
        edges = collect(0.0:0.5:100.0)
        for threads in (:off, :on)
            report!("fill", "$(EVENTS) events threads=$threads", measure(() -> histogram(x, edges; threads)))
        end
    end
    if "rootmap" in CASES
        h = Hist1D(collect(0.0:1.0:20.0), collect(1.0:20.0), collect(1.0:20.0))
        jobs = 1:8
        for workers in (1, 4)
            report!("rootmap", "8 PNG renders workers=$workers", measure(() ->
                rootmap(i -> render_root(h, GraphSpec(), joinpath(dir, "map-$workers-$i"); draw="E1"), jobs; workers)))
        end
    end
end

if !isempty(ARGS)
    open(ARGS[1], "w") do io
        println(io, "# ROOTPlots ", pkgversion(ROOTPlots), "\tJulia ", VERSION, "\tROOT ", built.version, "\t",
                Sys.KERNEL, "\t", Sys.ARCH, "\tthreads ", Threads.nthreads(:default))
        println(io, "case\tvariant\tmedian_ms\tiqr_ms\tn")
        for r in rows
            println(io, join((r[1], r[2], @sprintf("%.4f", r[3]), @sprintf("%.4f", r[4]), r[5]), '\t'))
        end
    end
end
