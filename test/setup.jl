using Test
using ROOTPlots

const ROOT_TEST_MODE = get(ENV, "ROOTPLOTS_TEST_ROOT", "auto")
ROOT_TEST_MODE in ("auto", "required", "off") ||
    error("ROOTPLOTS_TEST_ROOT must be auto, required, or off")
const TEST_NATIVE_ROOT = ROOT_TEST_MODE != "off" && root_available()
ROOT_TEST_MODE == "required" && !TEST_NATIVE_ROOT &&
    error("Native ROOT tests are required, but root-config was not found")
@info "ROOTPlots test mode" mode=ROOT_TEST_MODE native=TEST_NATIVE_ROOT

if !TEST_NATIVE_ROOT
    @testset "native ROOT integration (not run)" begin
        @test_skip false
    end
end

"""
    run_native(exe, args...; timeout=120)

Run a native backend executable directly (tests of its request validation), with
standard output discarded and standard error captured. The process is killed
after `timeout` seconds. Returns `(exitcode, termsignal, timed_out, stderr)`.
"""
function run_native(exe, args...; timeout=120)
    err = Pipe()
    process = run(pipeline(ignorestatus(Cmd([String(exe), map(string, args)...])); stdout=devnull, stderr=err);
                  wait=false)
    close(err.in)
    reader = @async read(err, String)
    timed_out = Ref(false)
    timer = Timer(timeout) do _
        if process_running(process)
            timed_out[] = true
            kill(process, Base.SIGKILL)
        end
    end
    wait(process)
    close(timer)
    return (exitcode=process.exitcode, termsignal=process.termsignal, timed_out=timed_out[],
            stderr=fetch(reader))
end

"""
    rejects_request(exe, request, stem, options...; timeout=120)

True when the native backend refuses `request` with a controlled validation
error: exit code 2, no signal, no timeout, and no output file for `stem`
(no `stem.png`, `stem.pdf`, or other `stem.*` file). `options` are further
command-line arguments (e.g. `"--formats=gif"`).
"""
function rejects_request(exe, request, stem, options...; timeout=120)
    r = run_native(exe, request, stem, options...; timeout)
    dir, base = dirname(stem), basename(stem)
    outputs = isdir(dir) ? filter(f -> startswith(f, base * "."), readdir(dir)) : String[]
    ok = r.exitcode == 2 && r.termsignal == 0 && !r.timed_out && isempty(outputs)
    ok || @warn "native request was not rejected cleanly" request r.exitcode r.termsignal r.timed_out outputs r.stderr
    return ok
end
