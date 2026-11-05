# Native process limits, cancellation, the root-config metadata cache (no ROOT
# needed: POSIX shell fixtures), and ROOT output formats / reproducibility records.
using ROOTPlots: run_backend, backend_timeout, timeout_setting, timeout_keyword,
                 BackendContext, BACKEND_CONTEXT, cancel!, root_config_metadata,
                 clear_backend_cache!, backend_identity, compile_root_backend
using Base.ScopedValues: with
import TOML

rt_sh(script) = Cmd(["/bin/sh", "-c", script])
rt_alive(pid) = ccall(:kill, Cint, (Cint, Cint), pid, 0) == 0
"Wait (bounded) until process pid no longer exists; true when it is gone."
function rt_gone(pid; within=10.0)
    deadline = time() + within
    while time() < deadline
        rt_alive(pid) || return true
        sleep(0.05)
    end
    return !rt_alive(pid)
end
rt_pid(path; within=10.0) = (timedwait(() -> isfile(path) && !isempty(strip(read(path, String))), within);
                              parse(Int, strip(read(path, String))))
rt_catch(f) = try f(); nothing catch e; e end
rt_open_fds() = Sys.islinux() ? length(readdir("/proc/self/fd")) : length(readdir("/dev/fd"))

@testset "backend time-limit settings (no ROOT)" begin
    withenv("ROOTPLOTS_BACKEND_TIMEOUT" => nothing) do
        @test backend_timeout() == 600.0
    end
    for (value, expected) in (("", 600.0), ("  ", 600.0), ("2.5", 2.5), ("30", 30.0),
                              ("0", Inf), ("0.0", Inf), ("none", Inf), ("NONE", Inf))
        withenv("ROOTPLOTS_BACKEND_TIMEOUT" => value) do
            @test backend_timeout() == expected
        end
    end
    for bad in ("abc", "-1", "inf", "nan", "10s", "1e400")
        withenv("ROOTPLOTS_BACKEND_TIMEOUT" => bad) do
            err = rt_catch(backend_timeout)
            @test err isa ArgumentError && occursin("ROOTPLOTS_BACKEND_TIMEOUT", err.msg)
        end
    end
    @test timeout_setting("ROOTPLOTS_ROOT_CONFIG_TIMEOUT_UNSET_FOR_TEST", 60.0) == 60.0
    withenv("ROOTPLOTS_ROOT_CONFIG_TIMEOUT" => nothing, "ROOTPLOTS_COMPILE_TIMEOUT" => nothing) do
        @test ROOTPlots.root_config_timeout() == 60.0 && ROOTPlots.compile_timeout() == 600.0
    end
    @test timeout_keyword(5) == 5.0 && timeout_keyword(Inf) == Inf
    for bad in (0, -1, NaN, "5", nothing)
        @test_throws ArgumentError timeout_keyword(bad)
    end
    @test_throws ArgumentError rootmap(identity, 1:2; timeout=0)
    @test_throws ArgumentError rootmap(identity, 1:2; timeout=-3)
    # Old six-field constructor still works (no time-limit event).
    e = RootBackendError("op", "", "exe", 1, 0, "")
    @test !e.timedout && e.timeout == Inf
end

if Sys.isunix()
@testset "backend runner: timeouts, process groups, interrupts (no ROOT)" begin
    mktempdir() do dir
        # Invalid settings fail before anything is started.
        marker = joinpath(dir, "started")
        withenv("ROOTPLOTS_BACKEND_TIMEOUT" => "soon") do
            @test_throws ArgumentError run_backend(rt_sh("touch '$marker'"); operation="guard")
            ran = Threads.Atomic{Int}(0)
            @test_throws ArgumentError rootmap(1:3) do i
                Threads.atomic_add!(ran, 1)
            end
            @test ran[] == 0
        end
        @test !isfile(marker)
        @test_throws ArgumentError run_backend(rt_sh("touch '$marker'"); operation="guard", timeout=0)
        @test !isfile(marker)

        # The child leads its own process group, never ours.
        out = run_backend(rt_sh("echo \$\$ \$(ps -o pgid= -p \$\$)"); operation="group", capture=true)
        pid, pgid = parse.(Int, split(out))
        @test pid == pgid
        @test pgid != ccall(:getpgrp, Cint, ())
        @test !ROOTPlots.owns_process_group(getpid())

        # A sleeping process is stopped (SIGTERM) after its limit; stderr is kept.
        elapsed = @elapsed err = rt_catch(() -> run_backend(rt_sh("echo 'sleeper started' >&2; sleep 120");
                                                            operation="sleep fixture", timeout=1.5))
        @test err isa RootBackendError && err.timedout && err.timeout == 1.5
        @test err.termsignal == 15 && err.executable == "sh"
        @test occursin("sleeper started", err.stderr)
        @test 1.5 <= elapsed < 30
        shown = sprint(showerror, err)
        @test occursin("native ROOT backend failed during sleep fixture: timed out after 1.5 s", shown)
        @test occursin("stopped by signal 15 (SIGTERM)", shown)

        # The whole owned group is terminated: no grandchild survives.
        pidfile = joinpath(dir, "grandchild.pid")
        err = rt_catch(() -> run_backend(rt_sh("sleep 300 & echo \$! > '$pidfile'; wait");
                                         operation="grandchild fixture", timeout=1.5))
        @test err isa RootBackendError && err.timedout
        grandchild = rt_pid(pidfile)
        @test rt_gone(grandchild)

        # A group ignoring SIGTERM gets SIGKILL after the 2 s grace.
        pidfile2 = joinpath(dir, "ignoring.pid")
        elapsed = @elapsed err = rt_catch(() -> run_backend(
            rt_sh("trap '' TERM; sleep 300 & echo \$! > '$pidfile2'; wait"); operation="stubborn", timeout=1.0))
        @test err isa RootBackendError && err.timedout && err.termsignal == 9
        @test occursin("SIGKILL", sprint(showerror, err))
        @test 1.0 + ROOTPlots.BACKEND_TERM_GRACE <= elapsed < 40
        @test rt_gone(rt_pid(pidfile2))

        # Large output while timing out: bounded head and tail, no deadlock.
        big = "head -c 3000000 /dev/zero | tr '\\0' 'x' >&2"
        err = rt_catch(() -> run_backend(rt_sh("echo FIRST >&2; $big; echo LAST >&2; sleep 120");
                                         operation="noisy", timeout=3.0))
        @test err isa RootBackendError && err.timedout
        @test startswith(err.stderr, "FIRST") && endswith(rstrip(err.stderr), "LAST")
        @test occursin("bytes of native output omitted", err.stderr) && sizeof(err.stderr) < 70_000

        # Disabled limits: 0 / none in the environment, Inf as a keyword.
        for value in ("0", "none")
            withenv("ROOTPLOTS_BACKEND_TIMEOUT" => value) do
                @test backend_timeout() == Inf
                @test run_backend(rt_sh("sleep 1; echo done"); operation="slow", capture=true) == "done\n"
            end
        end
        @test run_backend(rt_sh("echo fine"); operation="x", capture=true, timeout=Inf) == "fine\n"
        # A limit that is not reached changes nothing (exit code, stderr text).
        err = rt_catch(() -> run_backend(rt_sh("echo bad >&2; exit 4"); operation="x", timeout=60))
        @test err isa RootBackendError && !err.timedout && err.exitcode == 4 && err.timeout == 60.0
        @test !occursin("timed out", sprint(showerror, err))

        # Interrupting the waiting task stops the owned group and rethrows.
        pidfile3 = joinpath(dir, "interrupted.pid")
        task = @async run_backend(rt_sh("echo \$\$ > '$pidfile3.child'; sleep 300 & echo \$! > '$pidfile3'; wait");
                                  operation="interrupt fixture")
        child = rt_pid(pidfile3 * ".child"); grandchild = rt_pid(pidfile3)
        schedule(task, InterruptException(); error=true)
        result = rt_catch(() -> wait(task))
        @test result isa TaskFailedException && result.task.exception isa InterruptException
        @test rt_gone(child) && rt_gone(grandchild)

        # A descendant still holding the output streams after a successful exit is
        # stopped after the 2 s drain grace; the result is kept and a warning logged.
        pidfile4 = joinpath(dir, "holder.pid")
        elapsed = @elapsed out = @test_logs (:warn, r"left descendants holding its output streams") run_backend(
            rt_sh("sleep 300 >&2 & echo \$! > '$pidfile4'; echo result"); operation="holder", capture=true)
        @test out == "result\n" && elapsed < 30
        @test rt_gone(rt_pid(pidfile4))

        # An interrupt while the streams are drained (the child already exited, a
        # descendant still holds stderr) stops the owned group before the process
        # is released; the InterruptException is rethrown unchanged.
        pidfile5 = joinpath(dir, "drain.pid")
        task = @async run_backend(rt_sh("/bin/sleep 41 >&2 & echo \$! > '$pidfile5'; exit 0"); operation="drain")
        holder = rt_pid(pidfile5)
        sleep(0.7)                                  # inside the 2 s drain window
        @test !istaskdone(task) && rt_alive(holder)
        schedule(task, InterruptException(); error=true)
        result = rt_catch(() -> wait(task))
        @test result isa TaskFailedException && result.task.exception isa InterruptException
        @test rt_gone(holder; within=1.0)           # already stopped when run_backend returned
        sleep(1.0)
        @test !rt_alive(holder)

        # A failure of the stream reader is not reported as an interrupt, and the
        # owned group is still stopped.
        pidfile6 = joinpath(dir, "readfail.pid")
        h = ROOTPlots.launch_native(Cmd(rt_sh("sleep 300 & echo \$! > '$pidfile6'; wait"); detach=true), nothing)
        sleeper = rt_pid(pidfile6)
        ROOTPlots.stop_after_failure(h, ErrorException("reader failed"))
        @test h.reason == :failed && process_exited(h.process)
        @test rt_gone(sleeper)

        # No file descriptors leak across successful, failed and stopped runs.
        before = rt_open_fds()
        redirect_stderr(devnull) do
            for i in 1:10
                run_backend(rt_sh("echo out; echo err >&2"); operation="fd", capture=true)
                run_backend(rt_sh("echo err >&2"); operation="fd")
                rt_catch(() -> run_backend(rt_sh("exit 1"); operation="fd"))
            end
        end
        rt_catch(() -> run_backend(rt_sh("sleep 60"); operation="fd", timeout=0.5))
        GC.gc()
        @test rt_open_fds() <= before + 2
    end
end

@testset "rootmap: per-batch native limits and isolated cancellation (no ROOT)" begin
    mktempdir() do dir
        # One hung native job: indexed RootMapError with its diagnostic; the others succeed.
        done = Threads.Atomic{Int}(0)
        elapsed = @elapsed err = rt_catch() do
            rootmap(1:4; workers=4, timeout=2) do i
                if i == 2
                    run_backend(rt_sh("echo 'job 2 started' >&2; sleep 120"); operation="job 2")
                else
                    out = run_backend(rt_sh("echo $i"); operation="job $i", capture=true)
                    Threads.atomic_add!(done, 1)
                    parse(Int, out)
                end
            end
        end
        @test err isa RootMapError && first.(err.failures) == [2]
        inner = err.failures[1].second.ex
        @test inner isa RootBackendError && inner.timedout && inner.timeout == 2.0
        @test occursin("job 2 started", inner.stderr)
        @test occursin("[2] RootBackendError: native ROOT backend failed during job 2: timed out after 2 s",
                       sprint(showerror, err))
        @test done[] == 3 && elapsed < 40

        # The batch limit reaches worker tasks only; nested batches inherit or override it.
        seen = rootmap(1:2; workers=2, timeout=7) do i
            nested = rootmap(_ -> backend_timeout(), 1:2; workers=2)
            override = rootmap(_ -> backend_timeout(), 1:1; timeout=Inf)
            (backend_timeout(), nested, override)
        end
        @test all(s -> s == (7.0, [7.0, 7.0], [Inf]), seen)
        withenv("ROOTPLOTS_BACKEND_TIMEOUT" => "123") do
            @test backend_timeout() == 123.0
            @test rootmap(_ -> backend_timeout(), 1:2) == [123.0, 123.0]
        end

        # A cancelled scope starts nothing; a finished process leaves the registry.
        context = BackendContext(nothing, Inf)
        marker = joinpath(dir, "never")
        with(BACKEND_CONTEXT => context) do
            @test run_backend(rt_sh("true"); operation="x") === nothing
            @test isempty(context.processes)
            cancel!(context)
            @test_throws InterruptException run_backend(rt_sh("touch '$marker'"); operation="x")
        end
        @test !isfile(marker)

        # Interrupt batch A while batch B and an unrelated call run concurrently.
        started_a = Threads.Atomic{Int}(0)
        task_a = @async rootmap(1:6; workers=2) do i
            Threads.atomic_add!(started_a, 1)
            run_backend(rt_sh("echo \$\$ > '$dir/a-$i.pid'; sleep 300 & echo \$! > '$dir/a-$i.child'; wait");
                        operation="batch A job $i")
        end
        task_b = Threads.@spawn rootmap(1:3; workers=3) do i
            parse(Int, run_backend(rt_sh("echo \$\$ > '$dir/b-$i.pid'; sleep 4; echo $i"); operation="batch B", capture=true))
        end
        task_c = Threads.@spawn run_backend(rt_sh("sleep 4; echo unrelated"); operation="single", capture=true)
        @test timedwait(() -> all(i -> isfile("$dir/a-$i.child") && filesize("$dir/a-$i.child") > 0, 1:2) &&
                              all(i -> isfile("$dir/b-$i.pid"), 1:3), 30.0) == :ok
        a_pids = [rt_pid("$dir/a-$i.pid") for i in 1:2]
        a_children = [rt_pid("$dir/a-$i.child") for i in 1:2]
        # Interrupt the waiting batch task itself, as Ctrl-C would.
        elapsed = @elapsed begin
            schedule(task_a, InterruptException(); error=true)
            result = rt_catch(() -> wait(task_a))
        end
        @test result isa TaskFailedException && result.task.exception isa InterruptException
        @test elapsed < 30
        @test all(rt_gone, a_pids) && all(rt_gone, a_children)
        @test started_a[] == 2                          # no new job after the interrupt
        @test !any(i -> isfile("$dir/a-$i.pid"), 3:6)
        @test fetch(task_b) == [1, 2, 3]                # other batch unaffected
        @test fetch(task_c) == "unrelated\n"            # unrelated call unaffected

        # An InterruptException raised inside a job cancels that batch too.
        err = rt_catch() do
            rootmap(1:4; workers=2) do i
                i == 1 && throw(InterruptException())
                run_backend(rt_sh("sleep 300"); operation="sibling")
            end
        end
        @test err isa InterruptException

        # Nested: cancelling the outer batch stops the inner batch's processes.
        task_n = @async rootmap(1:1) do _
            rootmap(1:2; workers=2) do j
                run_backend(rt_sh("echo \$\$ > '$dir/n-$j.pid'; sleep 300"); operation="nested")
            end
        end
        @test timedwait(() -> all(j -> isfile("$dir/n-$j.pid") && filesize("$dir/n-$j.pid") > 0, 1:2), 30.0) == :ok
        nested_pids = [rt_pid("$dir/n-$j.pid") for j in 1:2]
        schedule(task_n, InterruptException(); error=true)
        result = rt_catch(() -> wait(task_n))
        @test result isa TaskFailedException && result.task.exception isa InterruptException
        @test all(rt_gone, nested_pids)
    end
end

@testset "root-config metadata cache and backend keys (no ROOT)" begin
    mktempdir() do dir
        counter = joinpath(dir, "calls")
        fake = joinpath(dir, "root-config")
        script(version) = "#!/bin/sh\necho \"\$*\" >> '$counter'\n" *
                          "case \"\$1\" in --version) echo $version;; *) echo '-I/fake/include -L/fake/lib -lCore';; esac\n"
        write(fake, script("6.99/01")); chmod(fake, 0o755)
        calls() = isfile(counter) ? countlines(counter) : 0
        clear_backend_cache!()
        m = root_config_metadata(fake)
        @test calls() == 2                                   # --version, --cflags --libs
        @test m.version == "6.99/01" && m.flags == ["-I/fake/include", "-L/fake/lib", "-lCore"]
        @test root_config_metadata(fake) == m && calls() == 2  # cache hit: nothing runs
        # Holding the build lock elsewhere does not block cached metadata.
        holder = Threads.@spawn lock(() -> sleep(3), ROOTPlots.ROOT_BUILD_LOCK)
        sleep(0.2)
        @test (@elapsed root_config_metadata(fake)) < 1.0
        wait(holder)
        # Another path to the same file, an edited file, and force all query again.
        alias = joinpath(dir, "alias-root-config")
        symlink(fake, alias)
        root_config_metadata(alias); @test calls() == 4
        root_config_metadata(alias); @test calls() == 4
        write(fake, script("6.100/02")); chmod(fake, 0o755)
        @test root_config_metadata(fake).version == "6.100/02" && calls() == 6
        root_config_metadata(fake; force=true); @test calls() == 8
        # Concurrent first lookups query once.
        clear_backend_cache!()
        rm(counter)
        rootmap(_ -> root_config_metadata(fake), 1:8; workers=8)
        @test calls() == 2

        # Backend keys: sources are hashed on every lookup, CXX resolved every call.
        source = joinpath(dir, "fixture.cpp"); header = joinpath(dir, "fixture.hpp")
        write(source, "int main() {}\n"); write(header, "// v1\n")
        compiler = joinpath(dir, "fake-cxx")
        write(compiler, "#!/bin/sh\nexit 0\n"); chmod(compiler, 0o755)
        compiler2 = joinpath(dir, "fake-cxx-2")
        cp(compiler, compiler2); chmod(compiler2, 0o755)
        withenv("CXX" => compiler) do
            n = calls()
            k1 = backend_identity(fake, source, "-std=c++17"; dependencies=[header]).key
            @test backend_identity(fake, source, "-std=c++17"; dependencies=[header]).key == k1
            write(header, "// v2\n")
            k2 = backend_identity(fake, source, "-std=c++17"; dependencies=[header]).key
            @test k2 != k1
            write(source, "int main() { return 0; }\n")
            k3 = backend_identity(fake, source, "-std=c++17"; dependencies=[header]).key
            @test k3 != k2
            @test backend_identity(fake, source, "-std=c++20"; dependencies=[header]).key != k3
            withenv("CXX" => compiler2) do
                @test backend_identity(fake, source, "-std=c++17"; dependencies=[header]).key != k3
            end
            withenv("CXX" => joinpath(dir, "missing-compiler")) do
                @test_throws ArgumentError backend_identity(fake, source, "-std=c++17")
            end
            @test calls() == n                               # keys never re-query root-config
        end

        # Compiler failures and time limits are RootBackendErrors; no staging is left.
        cache = joinpath(dir, "cache")
        failing = joinpath(dir, "failing-cxx")
        write(failing, "#!/bin/sh\necho 'fixture.cpp:1: error: intentional' >&2\nexit 1\n"); chmod(failing, 0o755)
        withenv("CXX" => failing) do
            err = rt_catch(() -> compile_root_backend(fake, source, "fixture", "-std=c++17"; build_dir=cache))
            @test err isa RootBackendError && err.exitcode == 1 && !err.timedout
            @test err.operation == "C++ compilation of fixture" && err.executable == "failing-cxx"
            @test occursin("error: intentional", err.stderr)
        end
        slow = joinpath(dir, "slow-cxx")
        write(slow, "#!/bin/sh\nsleep 120\n"); chmod(slow, 0o755)
        withenv("CXX" => slow, "ROOTPLOTS_COMPILE_TIMEOUT" => "1") do
            err = rt_catch(() -> compile_root_backend(fake, source, "fixture", "-std=c++17"; build_dir=cache))
            @test err isa RootBackendError && err.timedout && err.timeout == 1.0
        end
        cxx_marker = joinpath(dir, "compiler-ran")
        marking = joinpath(dir, "marking-cxx")
        write(marking, "#!/bin/sh\ntouch '$cxx_marker'\n"); chmod(marking, 0o755)
        withenv("CXX" => marking, "ROOTPLOTS_COMPILE_TIMEOUT" => "-5") do
            @test_throws ArgumentError compile_root_backend(fake, source, "fixture", "-std=c++17"; build_dir=cache)
        end
        @test !isfile(cxx_marker)
        @test all(!startswith(name, "compile-") for (_, dirs, _) in walkdir(cache) for name in dirs)
        # A hanging root-config query is bounded too.
        hanging = joinpath(dir, "hanging-root-config")
        write(hanging, "#!/bin/sh\nsleep 120\n"); chmod(hanging, 0o755)
        withenv("ROOTPLOTS_ROOT_CONFIG_TIMEOUT" => "1") do
            err = rt_catch(() -> root_config_metadata(hanging))
            @test err isa RootBackendError && err.timedout && err.operation == "root-config --version"
        end
        clear_backend_cache!()
    end
end
else
    @testset "backend runner fixtures need a POSIX shell" begin
        @test_skip false
    end
end

@testset "staged record files are removed when writing fails (no ROOT)" begin
    mktempdir() do dir
        destination = joinpath(dir, "out.provenance.toml")
        err = rt_catch(() -> ROOTPlots.staged_file(io -> (write(io, "partial"); error("injected writer failure")),
                                                  destination))
        @test err isa ErrorException && err.msg == "injected writer failure"   # not rewrapped
        @test isempty(readdir(dir))
        # A writer interrupted the same way is cleaned up and rethrown unchanged.
        err = rt_catch(() -> ROOTPlots.staged_file(io -> throw(InterruptException()), destination))
        @test err isa InterruptException && isempty(readdir(dir))
        # Success: one staged file with the written bytes, owner-only on POSIX.
        path = ROOTPlots.staged_file(io -> write(io, "ok"), destination)
        @test dirname(path) == dir && read(path, String) == "ok" && readdir(dir) == [basename(path)]
        Sys.isunix() && @test filemode(path) & 0o077 == 0
        rm(path)
    end
end

if TEST_NATIVE_ROOT
@testset "requested output formats only (native)" begin
    built = build_root()
    # ROOT's PDF info dictionary holds dates and the written file name (whose
    # length also shifts the cross-reference offsets): compare everything else.
    strip_pdf_metadata(s) = replace(first(split(s, "\nxref\n")), r"/(CreationDate|ModDate|Title) \([^)]*\)" => "")
    same_pdf(a, b) = strip_pdf_metadata(read(a, String)) == strip_pdf_metadata(read(b, String))
    mktempdir() do dir
        g = [GraphSeries([0., 1., 2.], [1., 3., 2.]; yerrors=[0.1, 0.2, 0.1], label="data")]
        h = Hist1D([0., 1., 2.], [2., 3.], [2., 3.])
        cmp = comparison(h, [Component("Model", h, "#336699")])
        requests = Dict(
            built.graphfit => render_root(g, GraphSpec(), joinpath(dir, "graph"); provenance=true).request,
            built.renderer => render_root([cmp], [PanelSpec("x")], joinpath(dir, "control"); provenance=true).request)
        for (exe, request) in requests
            name = basename(exe)
            outputs(stem) = sort(filter(f -> startswith(f, stem * "."), readdir(dir)))
            for (option, expected) in (("--formats=png", ["$name-png.png"]), ("--formats=pdf", ["$name-pdf.pdf"]),
                                       ("--formats=png,pdf", ["$name-both.pdf", "$name-both.png"]),
                                       ("--formats=pdf,png", ["$name-swapped.pdf", "$name-swapped.png"]))
                stem = "$name-" * Dict("--formats=png" => "png", "--formats=pdf" => "pdf",
                                       "--formats=png,pdf" => "both", "--formats=pdf,png" => "swapped")[option]
                r = run_native(exe, request, joinpath(dir, stem), option)
                @test r.exitcode == 0 && r.termsignal == 0 && !r.timed_out
                @test outputs(stem) == expected
            end
            # Legacy invocation without the option still writes both files.
            r = run_native(exe, request, joinpath(dir, "$name-legacy"))
            @test r.exitcode == 0 && outputs("$name-legacy") == ["$name-legacy.pdf", "$name-legacy.png"]
            # PNG alone is the same image as the PNG of a two-format run; the PDF
            # is equal apart from ROOT's dates and its file-name title.
            @test read(joinpath(dir, "$name-png.png")) == read(joinpath(dir, "$name-both.png"))
            @test read(joinpath(dir, "$name-png.png")) == read(joinpath(dir, "$name-legacy.png"))
            @test same_pdf(joinpath(dir, "$name-pdf.pdf"), joinpath(dir, "$name-both.pdf"))
            @test same_pdf(joinpath(dir, "$name-pdf.pdf"), joinpath(dir, "$name-legacy.pdf"))
            for bad in ("--formats=gif", "--formats=", "--format=png", "--formats=png,png", "--formats=PNG", "png")
                @test rejects_request(exe, request, joinpath(dir, "$name-bad"), bad)
            end
            @test rejects_request(exe, request, joinpath(dir, "$name-bad"), "--formats=png", "extra")
        end
        # Public API: keys unchanged, unrequested formats are nothing and absent.
        for (format, has_png, has_pdf) in ((:png, true, false), (:pdf, false, true), (:both, true, true))
            stem = joinpath(dir, "api-$format")
            p = render_root(g, GraphSpec(), stem; format)
            @test keys(p) == (:png, :pdf)
            @test (p.png !== nothing) == has_png == isfile(stem * ".png")
            @test (p.pdf !== nothing) == has_pdf == isfile(stem * ".pdf")
            c = render_root([cmp], [PanelSpec("x")], stem * "-control"; format)
            @test keys(c) == (:png, :pdf)
            @test isfile(stem * "-control.png") == has_png && isfile(stem * "-control.pdf") == has_pdf
        end
        @test read(joinpath(dir, "api-png.png")) == read(joinpath(dir, "api-both.png"))
        @test_throws ArgumentError render_root(g, GraphSpec(), joinpath(dir, "x"); format=:svg)

        # PDF metadata (no provenance): ROOT's /Title is the relative file name the
        # renderer saved inside its private directory, never a directory path.
        h2 = Hist2D([0., 1., 2.], [0., 1., 2.], [1. 3.; 2. 4.], fill(0.25, 2, 2))
        forbidden = filter(s -> length(s) >= 4, unique([dir, tempdir(), realpath(tempdir()), homedir(),
                                                         realpath(dirname(dir)), "/var/folders", "/tmp/"]))
        pdfs = [render_root(g, GraphSpec(), joinpath(dir, "meta-graph"); format=:pdf).pdf,
                render_root([cmp], [PanelSpec("x")], joinpath(dir, "meta-control"); format=:both).pdf,
                render_root(PlotFigure([PlotPanel(h, GraphSpec()), PlotPanel(h2, Graph2DSpec())]; layout=(1, 2)),
                            joinpath(dir, "meta-figure"); format=:pdf).pdf,
                render_root(h, GraphSpec(), joinpath(dir, "meta-hist"); draw="E1", format=:pdf).pdf]
        for pdf in pdfs
            text = read(pdf, String)
            titles = [m.captures[1] for m in eachmatch(r"/Title \(([^)]*)\)", text)]   # info + outlines
            @test "rendered.pdf" in titles && !any(t -> occursin('/', t), titles)
            @test !any(s -> occursin(s, text), forbidden)
        end
    end
end

@testset "opt-in reproducibility records (native)" begin
    built = build_root()
    mktempdir() do base
        # A private-looking directory name must not reach any record.
        marker = "rpPrivateMarker7Q"
        dir = joinpath(base, marker * "-outputs")
        mkpath(dir)
        label = "user label 42"                       # intentional request data
        h = Hist1D([0., 1., 2., 3., 4.], [4., 9., 7., 3.], [4., 9., 7., 3.])
        fit = fit_root(h, :pol1; initial=[5., 0.])
        h2 = Hist2D([0., 1., 2.], [0., 1., 2.], [1. 3.; 2. 4.], fill(0.25, 2, 2))
        g = [GraphSeries([0., 1., 2.], [1., 3., 2.]; yerrors=[0.1, 0.2, 0.1], label)]
        cmp = comparison(h, [Component("Model", h, "#336699")])
        panel = PlotPanel(h, GraphSpec(xlabel=label))
        nested = PlotFigure([panel, PlotFigure([panel, PlotPanel(h2, Graph2DSpec())]; layout=(2, 1))];
                            layout=(1, 2), widths=[1, 2])
        cases = [
            ("graph", "rootplots-graphfit", "DRAW", "CMS",
             stem -> render_root(g, GraphSpec(style=:CMS), stem; provenance=true, format=:both)),
            ("hist", "rootplots-graphfit", "DRAWHISTRAW", "",
             stem -> render_root(h, GraphSpec(), stem; draw="E1", provenance=true)),
            ("steps", "rootplots-graphfit", "DRAW", "",
             stem -> render_root(h, GraphSpec(), stem; provenance=true)),
            ("comparison", "rootplots-render", "SETTINGS", "",
             stem -> render_root([cmp], [PanelSpec(label)], stem; provenance=true, format=:both)),
            ("nested", "rootplots-graphfit", "DRAWLAYOUT", "",
             stem -> render_root(nested, stem; provenance=true)),
            ("fitpull", "rootplots-graphfit", "DRAWFIT", "",
             stem -> render_root(h, fit, GraphSpec(), stem; provenance=true)),
            ("th2", "rootplots-graphfit", "DRAWLAYOUT", "",
             stem -> render_root(h2, Graph2DSpec(), stem; draw="COLZ", provenance=true, format=:pdf)),
        ]
        forbidden = filter(s -> length(s) >= 4,
                           unique([marker, dir, base, tempdir(), homedir(), realpath(tempdir()),
                                   dirname(built.graphfit), ROOTPlots.root_cache_dir(),
                                   dirname(built.root_config), dirname(ROOTPlots.resolve_compiler())]))
        user = get(ENV, "USER", get(ENV, "USERNAME", ""))
        host = gethostname()
        for (name, executable, kind, style, render) in cases
            stem = joinpath(dir, name)
            p = render(stem)
            @test keys(p) == (:png, :pdf, :request, :provenance)
            @test p.request == stem * ".request.tsv" && p.provenance == stem * ".provenance.toml"
            @test isfile(p.request) && isfile(p.provenance)
            record = TOML.parsefile(p.provenance)
            @test record["record"] == "ROOTPlots provenance 1"
            @test record["rootplots"]["version"] == string(pkgversion(ROOTPlots))
            @test record["julia"]["version"] == string(VERSION)
            @test record["root"]["version"] == built.version
            @test !isempty(record["compiler"]["version"]) && !occursin('/', record["compiler"]["name"])
            @test record["system"]["kernel"] == string(Sys.KERNEL) && record["system"]["arch"] == string(Sys.ARCH)
            @test record["backend"]["executable"] == executable
            exe = executable == "rootplots-render" ? built.renderer : built.graphfit
            @test record["backend"]["key"] == basename(dirname(exe))
            @test record["request"]["file"] == basename(p.request)
            @test record["request"]["sha256"] == bytes2hex(open(ROOTPlots.sha256, p.request))
            @test record["request"]["style"] == style
            @test record["request"]["kind"] == kind
            formats = [f for (f, path) in (("png", p.png), ("pdf", p.pdf)) if path !== nothing]
            @test record["outputs"]["formats"] == formats
            for f in formats
                path = getfield(p, Symbol(f))
                @test record["outputs"][f]["file"] == basename(path)
                @test record["outputs"][f]["sha256"] == bytes2hex(open(ROOTPlots.sha256, path))
            end
            # Inline data only: no file references or directories in the request.
            request = read(p.request, String)
            @test !any(s -> occursin(s, request), forbidden)
            # The record holds file names only; no paths, user, host, or flags.
            text = read(p.provenance, String)
            @test !any(s -> occursin(s, text), forbidden)
            @test !occursin("-I", text) && !occursin("-L", text) && !occursin("\"/", text)
            isempty(user) || @test !occursin("/" * user, text) && !any(==(user), split(text, r"[\s\"=]+"))
            length(host) >= 3 && @test !occursin(host, text)
            # Replay: compiled backend + recorded request + chosen format.
            replay = joinpath(dir, name * "-replay")
            r = run_native(exe, p.request, replay, "--formats=png")
            @test r.exitcode == 0
            if p.png !== nothing
                @test read(replay * ".png") == read(p.png)
            else
                rendered_png = render(joinpath(dir, name * "-again")).request  # same request again
                @test read(rendered_png) == read(p.request)
                @test isfile(replay * ".png")
            end
        end
        # The graph label is user data and stays verbatim in the request.
        @test occursin(label, read(joinpath(dir, "graph.request.tsv"), String))
        @test occursin(label, read(joinpath(dir, "comparison.request.tsv"), String))
        @test occursin(label, read(joinpath(dir, "nested.request.tsv"), String))

        # Default: unchanged return value, no records.
        plain = render_root(g, GraphSpec(), joinpath(dir, "plain"))
        @test keys(plain) == (:png, :pdf) && plain.pdf === nothing
        @test sort(filter(f -> startswith(f, "plain."), readdir(dir))) == ["plain.png"]
        # A failed rendering leaves neither images nor records.
        withenv("ROOTPLOTS_BACKEND_TIMEOUT" => "0.01") do
            err = try render_root(g, GraphSpec(), joinpath(dir, "failed"); provenance=true); nothing catch e; e end
            @test err isa RootBackendError && err.timedout
        end
        @test isempty(filter(f -> startswith(f, "failed"), readdir(dir)))
        @test isempty(filter(f -> startswith(f, "jl_") || startswith(f, "tmp"), readdir(dir)))
        # A rootmap time limit applies to the native render of each job.
        err = try
            rootmap(1:2; timeout=0.01) do i
                render_root(g, GraphSpec(), joinpath(dir, "batch-$i"))
            end
            nothing
        catch e
            e
        end
        @test err isa RootMapError && first.(err.failures) == [1, 2]
        @test all(f -> f.second.ex isa RootBackendError && f.second.ex.timedout, err.failures)
        @test fetch(Threads.@spawn render_root(g, GraphSpec(), joinpath(dir, "after"))).png !== nothing
    end
end

@testset "relative TMPDIR (native, subprocess)" begin
    # A relative TMPDIR makes Julia's temporary directories relative; the renderer
    # runs in its private directory, so every path it gets must be absolute.
    built = build_root()
    mktempdir() do base
        mkpath(joinpath(base, "reltmp"))
        script = """
            using ROOTPlots
            g = [GraphSeries([0., 1., 2.], [1., 3., 2.]; yerrors=[0.1, 0.2, 0.1], label="data")]
            p = render_root(g, GraphSpec(), joinpath("out", "plot"); format=:both, provenance=true)
            println(join((p.png, p.pdf, p.request, p.provenance), '\\t'))
            println(ENV["TMPDIR"])
            """
        err = IOBuffer()
        cmd = setenv(`$(Base.julia_cmd()) --startup-file=no --project=$(pkgdir(ROOTPlots)) -e $script`,
                     merge(ENV, Dict("TMPDIR" => "reltmp")); dir=base)
        out = read(pipeline(ignorestatus(cmd); stderr=err), String)
        stderr_text = String(take!(err))
        lines = split(chomp(out), '\n')
        @test length(lines) == 2
        if length(lines) == 2
            png, pdf, request, record = split(lines[1], '\t')
            # Caller-relative output paths; Julia's ENV is not changed.
            @test (png, pdf, request, record) == ("out/plot.png", "out/plot.pdf", "out/plot.request.tsv",
                                                  "out/plot.provenance.toml")
            @test lines[2] == "reltmp"
            @test all(f -> isfile(joinpath(base, f)), (png, pdf, request, record))
            titles = [m.captures[1] for m in eachmatch(r"/Title \(([^)]*)\)", read(joinpath(base, pdf), String))]
            @test "rendered.pdf" in titles && !any(t -> occursin('/', t), titles)
            @test !occursin(base, read(joinpath(base, pdf), String))
            parsed = TOML.parsefile(joinpath(base, record))
            @test parsed["outputs"]["png"]["sha256"] == bytes2hex(open(ROOTPlots.sha256, joinpath(base, png)))
            r = run_native(built.graphfit, joinpath(base, request), joinpath(base, "replayed"), "--formats=png")
            @test r.exitcode == 0 && read(joinpath(base, "replayed.png")) == read(joinpath(base, png))
        end
        @test !occursin("temporary", lowercase(stderr_text)) && !occursin("Error", stderr_text)
        @test isempty(readdir(joinpath(base, "reltmp")))        # private directories removed
    end
end
end
