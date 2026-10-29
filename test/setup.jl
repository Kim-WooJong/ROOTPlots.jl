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
