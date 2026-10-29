@testset "local ROOT compilation and cache reuse" begin
    withenv("ROOT_CONFIG" => nothing, "ROOTSYS" => nothing, "PATH" => "") do
        @test_throws ArgumentError build_root()
    end
    if TEST_NATIVE_ROOT
        # Exercise real compilation/linking rather than only finding root-config.
        built = build_root(; force=true)
        @test built.version == strip(read(Cmd([built.root_config, "--version"]), String))
        @test isfile(built.renderer) && isfile(built.graphfit)
        stamps = (mtime(built.renderer), mtime(built.graphfit))
        @test build_root(; root_config=built.root_config) == built
        @test (mtime(built.renderer), mtime(built.graphfit)) == stamps
        mktempdir() do dir
            data = GraphSeries([0., 1., 2., 3.], [1., 3., 5., 7.];
                               yerrors=fill(0.2, 4))
            fit = fit_root(data, :pol1; initial=[0., 0.], root_config=built.root_config)
            @test fit.status == 0 && fit.valid
            @test fit.parameters ≈ [1., 2.] atol=1e-6
            files = render_root([fit_curve(fit), data], GraphSpec(),
                                joinpath(dir, "compiled"); format=:both,
                                root_config=built.root_config)
            @test read(files.png)[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
            @test startswith(read(files.pdf, String), "%PDF-")
            h = Hist1D([0., 1., 2.], [2., 3.], [2., 3.])
            cmp = comparison(h, [Component("Model", h, "#336699")])
            control = render_root([cmp], [PanelSpec("x")], joinpath(dir, "control");
                                  format=:both, root_config=built.root_config)
            @test isfile(control.png) && isfile(control.pdf)
        end
    end
end
