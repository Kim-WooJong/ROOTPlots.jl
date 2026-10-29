using ROOTPlots

println("ROOTPlots environment check")
println("Julia: ", VERSION)
println("Compiler: ", something(Sys.which(get(ENV, "CXX", "c++")), "NOT FOUND"))
root_available() || error("root-config not found; set ROOT_CONFIG to its executable path or add ROOT/bin to PATH")
built = build_root()
println("ROOT: ", built.version, " (", built.root_config, ")")
println("Control backend: ", built.renderer)
println("Graph/fit backend: ", built.graphfit)

# Build success alone does not prove that the dynamic loader or rendering works.
mktempdir() do dir
    data = GraphSeries([0., 1., 2.], [1., 3., 5.]; yerrors=fill(.2, 3))
    fit = fit_root(data, :pol1; initial=[0., 1.])
    fit.status == 0 && fit.valid && isapprox(fit.parameters, [1., 2.]; atol=1e-6) ||
        error("native fit smoke check failed")
    files = render_root([fit_curve(fit), data], GraphSpec(), joinpath(dir, "check"); format=:both)
    read(files.png)[1:8] == UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a] || error("invalid PNG")
    startswith(read(files.pdf, String), "%PDF-") || error("invalid PDF")
end
println("PASS: compile, link, native fit, PNG and PDF output")
