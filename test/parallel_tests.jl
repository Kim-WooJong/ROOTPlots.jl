@testset "bounded threaded rootmap ordering and failure lifetime" begin
    @test rootmap(identity, Int[]) == []
    @test_throws ArgumentError rootmap(identity, 1:3; workers=0)
    @test_throws ArgumentError rootmap(identity, Int[]; workers=-1)
    state_lock = ReentrantLock()
    active = Ref(0)
    peak = Ref(0)
    thread_ids = Set{Int}()
    results = rootmap(1:9; workers=3) do i
        lock(state_lock) do
            active[] += 1
            peak[] = max(peak[], active[])
            push!(thread_ids, Threads.threadid())
        end
        try
            sleep(0.005 * (10-i)) # deliberately reverse completion pressure
            i^2
        finally
            lock(state_lock) do
                active[] -= 1
            end
        end
    end
    @test results == collect(1:9).^2
    @test peak[] == 3
    @test active[] == 0
    @info "rootmap worker threads observed" default_threads=Threads.nthreads(:default) thread_ids
    @test rootmap(string, 1:3; workers=1) == ["1", "2", "3"]
    @test rootmap(identity, (i for i in 1:2); workers=100) == [1, 2]
    completed = Threads.Atomic{Int}(0)
    error = try
        rootmap(1:5; workers=3) do i
            try
                sleep(0.005i)
                i in (2, 4) && throw(ArgumentError("bad job $i"))
                i
            finally
                Threads.atomic_add!(completed, 1)
            end
        end
        nothing
    catch exception
        exception
    end
    @test error isa RootMapError
    @test completed[] == 5 # no detached workers after error
    @test first.(error.failures) == [2, 4]
    @test all(p -> p.second.ex isa ArgumentError, error.failures)
    @test occursin("bad job 4", sprint(showerror, error))
end

@testset "concurrent native ROOT fits and outputs match serial inputs" begin
    if TEST_NATIVE_ROOT
        x = collect(0.:0.5:5.)
        datasets = [GraphSeries(x, a .+ (a+1).*x; yerrors=fill(0.2, length(x))) for a in 1:4]
        serial = [fit_root(data, :pol1; initial=[0., 0.]) for data in datasets]
        parallel = rootmap(datasets; workers=4) do data
            fit_root(data, :pol1; initial=[0., 0.])
        end
        @test all(f -> f.status == 0 && f.valid, parallel)
        @test all(i -> parallel[i].parameters ≈ serial[i].parameters, eachindex(serial))
        @test all(i -> parallel[i].covariance ≈ serial[i].covariance, eachindex(serial))
        mktempdir() do dir
            outputs = rootmap(1:4; workers=4) do i
                stem = joinpath(dir, "parallel-$i")
                if i == 1
                    h = Hist1D([0., 1., 2.], [2., 3.], [2., 3.])
                    cmp = comparison(h, [Component("Model", h, "#336699")])
                    render_root([cmp], [PanelSpec("x")], stem; format=:both)
                elseif i == 2
                    grid = Grid2D([0., 1., 2.], [0., 1., 2.], [1. 2.; 3. 4.])
                    render_root(grid, Graph2DSpec(), stem; format=:both)
                elseif i == 3
                    h = Hist1D([0., 1., 2., 3.], [2., 4., 3.], [2., 4., 3.])
                    fit = fit_root(h, :pol0; initial=[1.])
                    fit.status == 0 && fit.valid || error("parallel histogram fit failed")
                    render_root(h, fit, GraphSpec(), stem; format=:both)
                else
                    render_root([fit_curve(parallel[i]), datasets[i]], GraphSpec(), stem; format=:both)
                end
            end
            @test length(unique(p.png for p in outputs)) == 4
            @test all(p -> read(p.png)[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a], outputs)
            @test all(p -> startswith(read(p.pdf, String), "%PDF-"), outputs)
            # Identical graph inputs rendered sequentially/parallel have identical
            # PNG bytes under this ROOT installation (PDF metadata is not compared).
            reference = render_root([fit_curve(parallel[4]), datasets[4]], GraphSpec(),
                                    joinpath(dir, "serial-reference"))
            @test read(reference.png) == read(outputs[4].png)
        end
    end
end

@testset "cold/forced builds publish complete configuration-keyed artifacts" begin
    if TEST_NATIVE_ROOT
        config = ROOTPlots.find_root_config()
        mktempdir() do dir
            source = joinpath(dir, "minimal.cpp")
            header = joinpath(dir, "value.h")
            write(header, "#define OUTPUT 1\n")
            write(source, "#include <iostream>\n#include \"value.h\"\nint main() { std::cout << OUTPUT; }\n")
            cache = joinpath(dir, "cache")
            build(; force=false) = ROOTPlots.compile_root_backend(config, source,
                "fixture", "-std=c++17"; force, build_dir=cache)
            paths = rootmap(_ -> build(), 1:4; workers=4)
            @test length(unique(paths)) == 1
            first_exe = first(paths)
            @test read(Cmd([first_exe]), String) == "1"
            stamp = mtime(first_exe)
            @test build() == first_exe && mtime(first_exe) == stamp
            @test rootmap(1:3; workers=3) do _
                read(Cmd([build(; force=true)]), String)
            end == fill("1", 3)
            # Header changes are not in the source key: force must preserve the
            # exact same cached executable if the compiler fails.
            write(header, "#error intentional force build failure\n")
            redirect_stderr(devnull) do
                @test_throws ProcessFailedException build(; force=true)
            end
            @test read(Cmd([first_exe]), String) == "1"
            write(source, "#include <iostream>\nint main() { std::cout << 2; }\n")
            second_exe = build()
            @test second_exe != first_exe
            @test read(Cmd([first_exe]), String) == "1" # old key remains intact
            @test read(Cmd([second_exe]), String) == "2"
            alternate = ROOTPlots.compile_root_backend(config, source,
                "fixture", "-std=c++20"; build_dir=cache)
            @test alternate != second_exe
            @test all(!startswith(name, "compile-") for (_, dirs, _) in walkdir(cache) for name in dirs)
            # Separate Julia processes do not share ROOT_BUILD_LOCK. Atomic
            # publication must also tolerate their simultaneous cold/force builds.
            crosscache = joinpath(dir, "cross-process")
            script = "using ROOTPlots; exe = ROOTPlots.compile_root_backend(" *
                repr(config) * ", " * repr(source) * ", \"fixture\", \"-std=c++17\"; " *
                "build_dir=" * repr(crosscache) * ", force=true); print(read(Cmd([exe]), String))"
            command = `$(Base.julia_cmd()) --threads=1 --startup-file=no --project=$(dirname(Base.active_project())) -e $script`
            processes = [run(pipeline(command; stdout=joinpath(dir,"child-$i.out"),
                                      stderr=joinpath(dir,"child-$i.err")); wait=false) for i in 1:2]
            foreach(wait, processes)
            @test all(success, processes)
            @test all(i -> read(joinpath(dir,"child-$i.out"), String) == "2", 1:2)
        end
    end
end
