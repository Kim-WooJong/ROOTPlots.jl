# Separate integration probe: verifies package-source writes are unnecessary.
using Pkg

mktempdir() do temporary
    checkout = joinpath(temporary, "readonly-source")
    mkpath(checkout)
    for entry in ("src", "backend", "Project.toml")
        cp(joinpath(@__DIR__, "..", entry), joinpath(checkout, entry))
    end
    for (root, dirs, files) in walkdir(checkout)
        foreach(file -> chmod(joinpath(root, file), 0o444), files)
        foreach(dir -> chmod(joinpath(root, dir), 0o555), dirs)
    end
    chmod(checkout, 0o555)
    try
        environment = joinpath(temporary, "analysis")
        Pkg.activate(environment)
        Pkg.develop(path=checkout)
        Pkg.instantiate()
        # New process mirrors real installation and avoids loading newly created
        # package bindings inside an already compiled closure (Julia world ages).
        probe = raw"""
        using Pkg, Test, ROOTPlots
        checkout = dirname(dirname(pathof(ROOTPlots)))
        @testset "clean analysis environment with read-only package source" begin
            names = Set(info.name for info in values(Pkg.dependencies()))
            @test isempty(intersect(names, Set(["CairoMakie", "Makie", "Plots"])))
            @test issubset(Set(["CSV", "DataFrames"]), names)
            @test !ispath(joinpath(checkout, "build"))
            built = build_root()
            @test isfile(built.renderer) && isfile(built.graphfit)
            @test !startswith(built.graphfit, checkout)
            @test !ispath(joinpath(checkout, "Manifest.toml"))
            data = GraphSeries([0., 1., 2.], [1., 3., 5.]; yerrors=fill(.2, 3))
            fit = fit_root(data, :pol1; initial=[0., 1.])
            @test fit.status == 0 && fit.valid
            @test fit.parameters ≈ [1., 2.] atol=1e-6
            files = render_root([data], GraphSpec(), joinpath(dirname(checkout), "plot"); format=:both)
            @test read(files.png)[1:8] == UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a]
            @test startswith(read(files.pdf, String), "%PDF-")
            @test !ispath(joinpath(checkout, "build"))
        end
        """
        run(`$(Base.julia_cmd()) --startup-file=no --project=$environment -e $probe`)
    finally
        # Restore permissions only on the disposable fixture, for cleanup.
        for (root, _, files) in walkdir(checkout)
            chmod(root, 0o755)
            foreach(file -> chmod(joinpath(root, file), 0o644), files)
        end
    end
end
