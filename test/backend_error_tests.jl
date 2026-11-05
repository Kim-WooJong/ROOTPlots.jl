# Native process failures keep their diagnostics (RootBackendError), success forwards stderr.
using ROOTPlots: run_backend

@testset "backend runner: native failure text, exit code, signal (no ROOT)" begin
    if Sys.isunix()
        sh(script) = Cmd(["/bin/sh", "-c", script])
        err = try
            run_backend(sh("echo 'Error in <TH1D::Fit>: bad input' >&2; exit 3"); operation="test operation")
            nothing
        catch e
            e
        end
        @test err isa RootBackendError
        @test err.exitcode == 3 && err.termsignal == 0 && err.operation == "test operation"
        @test err.executable == "sh" && err.request == ""
        @test occursin("Error in <TH1D::Fit>: bad input", err.stderr)
        shown = sprint(showerror, err)
        @test occursin("native ROOT backend failed during test operation: exit code 3", shown)
        @test occursin("\n    Error in <TH1D::Fit>: bad input", shown)
        # A crash (signal) is reported as such; text written before it is kept.
        crash = try run_backend(sh("echo 'about to crash' >&2; kill -KILL \$\$"); operation="crash"); nothing catch e; e end
        @test crash isa RootBackendError && crash.termsignal == 9
        @test occursin("about to crash", crash.stderr)
        @test occursin("terminated by signal 9 (SIGKILL)", sprint(showerror, crash))
        # Request context: first record of the request file; non-ASCII requests never run.
        mktempdir() do dir
            request = joinpath(dir, "request.tsv")
            write(request, "ROOTPlotsGraph\t1\nDRAWFIT\t1\n")
            e = try run_backend(sh("exit 2"); operation="fit/pull rendering", request); nothing catch x; x end
            @test e.request == "DRAWFIT" && occursin("(request DRAWFIT)", sprint(showerror, e))
            @test occursin("native stderr: (empty)", sprint(showerror, e))
            marker = joinpath(dir, "ran")
            write(request, "ROOTPlotsGraph\t1\nLABELS\tμ\n")
            @test_throws ArgumentError run_backend(sh("touch '$marker'"); operation="guard", request)
            @test !isfile(marker)
        end
        # Large output cannot block the child; failures keep a bounded head and tail.
        big = "head -c 3000000 /dev/zero | tr '\\0' 'x' >&2"
        big_err = try run_backend(sh("echo FIRST >&2; $big; echo LAST >&2; exit 1"); operation="big"); nothing catch e; e end
        @test big_err isa RootBackendError
        @test startswith(big_err.stderr, "FIRST") && endswith(rstrip(big_err.stderr), "LAST")
        @test occursin("bytes of native output omitted", big_err.stderr) && sizeof(big_err.stderr) < 70_000
        # Success: stdout capture works and stderr (ROOT warnings) is forwarded to Julia's stderr.
        mktempdir() do dir
            log = joinpath(dir, "stderr.txt")
            out = open(log, "w") do io
                redirect_stderr(io) do
                    run_backend(sh("echo 'Warning in <TFile>: kept' >&2; $big; echo RESULT"); operation="ok", capture=true)
                end
            end
            @test out == "RESULT\n"
            text = read(log, String)
            @test startswith(text, "Warning in <TFile>: kept") && sizeof(text) > 3_000_000
        end
        # rootmap keeps every job's own native diagnostic.
        e = try
            rootmap(1:4; workers=2) do i
                iseven(i) ? run_backend(sh("echo 'job $i failed natively' >&2; exit $i"); operation="job $i") : i
            end
            nothing
        catch x
            x
        end
        @test e isa RootMapError && first.(e.failures) == [2, 4]
        @test all(p -> p.second.ex isa RootBackendError, e.failures)
        @test [p.second.ex.exitcode for p in e.failures] == [2, 4]
        @test occursin("job 4 failed natively", e.failures[2].second.ex.stderr)
        shown = sprint(showerror, e)
        @test occursin("[2] RootBackendError: native ROOT backend failed during job 2: exit code 2", shown)
        @test occursin("\n          job 4 failed natively", shown)   # indented under its index
    else
        @test_skip false
    end
end

@testset "fit_root size limits match the native backend (no ROOT)" begin
    n = 1_000_001
    big = Hist1D(collect(0.0:n), ones(n), ones(n))
    err = try fit_root(big, :pol0; initial=[1.0]); nothing catch e; e end
    @test err isa ArgumentError && occursin("1,000,000 histogram bins", err.msg)
    g = GraphSeries(collect(1.0:n), ones(n); yerrors=ones(n))
    err = try fit_root(g, :pol0; initial=[1.0]); nothing catch e; e end
    @test err isa ArgumentError && occursin("1,000,000 graph points", err.msg)
end

if TEST_NATIVE_ROOT
    @testset "native validation failures arrive as RootBackendError" begin
        overflow = Profile1D([0.0, 1.0], [0.2, 0.4], [1e154, 1e154]; error=:gaussian)
        err = try profile_stats(overflow); nothing catch e; e end
        @test err isa RootBackendError
        @test err.operation == "profile statistics" && err.request == "PROFILESTATS"
        @test err.exitcode == 2 && err.executable == "rootplots-graphfit"
        @test occursin("ROOTPlots graphfit: profile aggregate statistics overflow", err.stderr)
        fine = Profile1D([0.0, 1.0], [0.2, 0.4], [1.0, 2.0])
        e = try rootmap(p -> profile_stats(p).mean, [fine, overflow, fine]; workers=3); nothing catch x; x end
        @test e isa RootMapError && first.(e.failures) == [2]
        @test e.failures[1].second.ex isa RootBackendError
        @test occursin("profile aggregate statistics overflow", sprint(showerror, e))
        # ROOT's messages on success still reach Julia's stderr (ROOTPLOTS_ROOT_VERBOSE=1).
        mktempdir() do dir
            log = joinpath(dir, "stderr.txt")
            withenv("ROOTPLOTS_ROOT_VERBOSE" => "1") do
                open(log, "w") do io
                    redirect_stderr(io) do
                        render_root([GraphSeries([0.0, 1.0], [1.0, 2.0])], GraphSpec(), joinpath(dir, "verbose"))
                    end
                end
            end
            @test occursin("TCanvas::Print", read(log, String))
        end
    end
end
