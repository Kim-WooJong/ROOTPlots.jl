@testset "DataFrames adapter and ROOT selection" begin
    data = DataFrame(x=[.2, 1.2], y=[1.2, .2], w=[2., -3.])
    edges = [0., 1., 2.]
    @test histogram(view(data, :, :), :x, edges; weight=:w).hist.sumw == [2., -3.]
    @test histogram(data, (:x, :y), (edges, edges); weight=:w).hist.sumw == [0. 2.; -3. 0.]
    @test_throws ArgumentError histogram(data, :absent, edges)
    @test_throws ArgumentError histogram(data, :x, edges; weight=:absent)
    @test_throws ArgumentError histogram(data, (:x, :absent), (edges, edges))
    @test ROOTPlots.find_root_config("/nonexistent/root-config") === nothing
    withenv("ROOT_CONFIG" => "/nonexistent/root-config") do
        @test !root_available()
        @test_throws ArgumentError build_root()
    end
    if TEST_NATIVE_ROOT
        cache = ROOTPlots.root_cache_dir()
        @test !startswith(realpath(cache), realpath(dirname(@__DIR__)) * string(Base.Filesystem.path_separator))
        built = build_root()
        @test startswith(built.graphfit, cache)
        withenv("CXX" => "/nonexistent/compiler") do
            @test_throws ArgumentError build_root()
        end
    end
end
