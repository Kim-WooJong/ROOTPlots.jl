# Native request validation: malformed requests end in exit code 2 without a signal
# or an image (helpers run_native/rejects_request in setup.jl). Covers the
# comparison renderer (rootplots-render) and the CMSCANVAS rows of rootplots-graphfit.

"""Replace the first match of a request line pattern, failing loudly when absent."""
function mutate(valid, pattern, replacement)
    occursin(pattern, valid) || error("test request has no match for $pattern")
    return replace(valid, pattern => replacement; count=1)
end

"""The request with every line whose tag is in tags removed."""
drop_rows(valid, tags...) = join(filter(l -> !(first(split(l, '\t')) in tags), split(valid, '\n'; keepempty=true)), '\n')

"""A tab-separated line with some (1-based) fields replaced."""
function set_fields(line, pairs...)
    f = split(line, '\t')
    for (k, v) in pairs
        f[k] = v
    end
    return join(f, '\t')
end

"""The line of the request that starts with prefix."""
request_line(valid, prefix) = first(filter(l -> startswith(l, prefix), split(valid, '\n')))

@testset "native validation: comparison renderer" begin
    if TEST_NATIVE_ROOT
        data = Hist1D([0.0, 1, 2, 3], [4.0, 6, 5], [4.0, 6, 5])
        mc = Hist1D([0.0, 1, 2, 3], [3.5, 6.5, 4.0], [0.3, 0.4, 0.3])
        cmp = comparison(data, [Component("Signal", 0.5 * mc, "#cc3311"), Component("Background", 0.5 * mc, "#4477aa")])
        spec = PanelSpec("m [GeV]"; lumi="1 fb^{-1}", style=:CMS,
                         xaxis=AxisStyle(title_size=0.05), canvas=CanvasStyle(grid_y=true))
        exe = ROOTPlots.build_root_backend(ROOTPlots.find_root_config())
        mktempdir() do dir
            input = joinpath(dir, "figure.tsv")
            ROOTPlots.write_root_spec(input, [cmp, cmp], [spec, spec], (0.5, 1.5))
            valid = read(input, String)
            r = run_native(exe, input, joinpath(dir, "ok"))
            @test r.exitcode == 0 && r.termsignal == 0 && isfile(joinpath(dir, "ok.png"))
            bin1 = request_line(valid, "BIN\t1\t1\t")
            bin2 = request_line(valid, "BIN\t1\t2\t")
            comp1 = request_line(valid, "COMP\t1\t1\t")
            axis1 = request_line(valid, "AXIS\t1")
            f = split(bin1, '\t')
            with_field(k, v) = join([f[1:k-1]; v; f[k+1:end]], '\t')
            cases = [
                "no panels" => drop_rows(valid, "PANEL", "AXIS", "EXPLABEL", "PCANVAS", "PAXIS", "BIN", "COMP"),
                "unknown row" => mutate(valid, "SETTINGS\t", "SETTING\t"),
                "unknown trailing row" => valid * "EXTRA\t1\n",
                "blank row" => mutate(valid, "\nPANEL\t2", "\n\nPANEL\t2"),
                "bad schema" => mutate(valid, "ROOTPlots\t2", "ROOTPlots\t3"),
                "duplicate SETTINGS" => mutate(valid, "STYLE\t", request_line(valid, "SETTINGS") * "\nSTYLE\t"),
                "duplicate STYLE" => mutate(valid, "ROOTSTYLE\t", request_line(valid, "STYLE\t") * "\nROOTSTYLE\t"),
                "duplicate ROOTSTYLE" => mutate(valid, "ROOTSTYLE\tCMS", "ROOTSTYLE\tCMS\nROOTSTYLE\tCMS"),
                "duplicate AXIS" => mutate(valid, "\nEXPLABEL", "\n" * request_line(valid, "AXIS\t1") * "\nEXPLABEL"),
                "duplicate EXPLABEL" => mutate(valid, "\nPCANVAS\t1", "\n" * request_line(valid, "EXPLABEL") * "\nPCANVAS\t1"),
                "duplicate PCANVAS" => mutate(valid, "\nPAXIS\t1\tx", "\n" * request_line(valid, "PCANVAS\t1") * "\nPAXIS\t1\tx"),
                "duplicate PAXIS" => mutate(valid, "\nPAXIS\t1\ty", "\n" * request_line(valid, "PAXIS\t1\tx") * "\nPAXIS\t1\ty"),
                "duplicate BIN" => mutate(valid, bin1, bin1 * "\n" * bin1),
                "duplicate PANEL id" => mutate(valid, "PANEL\t2\t", "PANEL\t1\t"),
                "panel id out of order" => mutate(valid, "PANEL\t1\t", "PANEL\t2\t"),
                "BIN index skips" => mutate(valid, bin2, replace(bin2, "BIN\t1\t2\t" => "BIN\t1\t3\t")),
                "BIN panel id" => mutate(valid, bin1, replace(bin1, "BIN\t1\t1\t" => "BIN\t2\t1\t")),
                "missing BIN" => mutate(valid, bin2 * "\n", ""),
                "COMP before BIN" => mutate(mutate(valid, comp1 * "\n", ""), bin1, comp1 * "\n" * bin1),
                "BIN after COMP" => mutate(mutate(valid, bin2 * "\n", ""), comp1, comp1 * "\n" * bin2),
                "COMP index" => mutate(valid, "COMP\t1\t1\t", "COMP\t1\t2\t"),
                "missing COMP" => drop_rows(valid, "COMP"),
                "negative component" => mutate(valid, comp1, join([split(comp1, '\t')[1:end-1]; "-1.0"], '\t')),
                "component color" => mutate(valid, "\t#cc3311\t", "\tred\t"),
                "empty component label" => mutate(valid, "\tSignal\t", "\t\t"),
                "COMP field count" => mutate(valid, comp1, comp1 * "\t1.0"),
                "BIN field count" => mutate(valid, bin1, join(f[1:end-1], '\t')),
                "BIN extra field" => mutate(valid, bin1, bin1 * "\t0"),
                "nonfinite content" => mutate(valid, bin1, with_field(6, "nan")),
                "infinite edge" => mutate(valid, bin1, with_field(4, "-inf")),
                "negative data error" => mutate(valid, bin1, with_field(7, "-0.5")),
                "negative prediction" => mutate(valid, bin1, with_field(9, "-1.0")),
                "negative ratio band" => mutate(valid, bin1, with_field(16, "-0.1")),
                "validity flag" => mutate(valid, bin1, with_field(17, "2")),
                "reversed edges" => mutate(valid, bin1, with_field(5, "-1.0")),
                "ratio range" => mutate(valid, request_line(valid, "SETTINGS"), "SETTINGS\t1.5\t0.5"),
                "nonfinite ratio range" => mutate(valid, request_line(valid, "SETTINGS"), "SETTINGS\t0.5\tinf"),
                "marker style" => mutate(valid, request_line(valid, "STYLE\t"), replace(request_line(valid, "STYLE\t"), r"^STYLE\t\d+" => "STYLE\t99")),
                "AXIS limits" => mutate(valid, axis1, set_fields(axis1, 6 => "10.0", 7 => "1.0")),
                "AXIS boolean" => mutate(valid, "AXIS\t1\t0\t", "AXIS\t1\t7\t"),
                "PCANVAS margin option" => mutate(valid, request_line(valid, "PCANVAS\t1"), replace(request_line(valid, "PCANVAS\t1"), r"^PCANVAS\t1\t-" => "PCANVAS\t1\t0.2")),
                "unknown style" => mutate(valid, "ROOTSTYLE\tCMS", "ROOTSTYLE\tNoSuchStyle"),
            ]
            for (name, broken) in cases
                @test broken != valid
                write(input, broken)
                ok = rejects_request(exe, input, joinpath(dir, "bad"))
                ok || @info "comparison case not rejected" name
                @test ok
            end
        end
    else
        @test_skip false
    end
end

@testset "native validation: CMSCANVAS rows" begin
    if TEST_NATIVE_ROOT
        exe = ROOTPlots.graphfit_executable(ROOTPlots.find_root_config())
        s = GraphSeries([1.0, 2.0, 3.0], [2.0, 3.0, 2.5])
        mktempdir() do dir
            input = joinpath(dir, "graph.tsv")
            spec = GraphSpec(style=:CMS, cms_label=(extra_space=0.05, y_title_offset=1.1))
            ROOTPlots.write_graph_spec(input, [s], spec, ROOTPlots.graph_limits([s], spec))
            valid = read(input, String)
            @test occursin("ROOTSTYLE\tCMS\nCMSCANVAS\t0.05\t1.1\n", valid)
            r = run_native(exe, input, joinpath(dir, "ok"))
            @test r.exitcode == 0 && r.termsignal == 0 && isfile(joinpath(dir, "ok.png"))
            atlas = GraphSpec(style=:ATLAS)
            ROOTPlots.write_graph_spec(input, [s], atlas, ROOTPlots.graph_limits([s], atlas))
            atlas_request = read(input, String)
            unstyled = mutate(atlas_request, "ROOTSTYLE\tATLAS\n", "")   # a scene without a style
            cases = [
                "stray (no style)" => mutate(unstyled, "SCENE\t1\n", "SCENE\t1\nCMSCANVAS\t0\t-\n"),
                "stray (ATLAS style)" => mutate(atlas_request, "ROOTSTYLE\tATLAS\n", "ROOTSTYLE\tATLAS\nCMSCANVAS\t0\t-\n"),
                "without its ROOTSTYLE" => mutate(valid, "ROOTSTYLE\tCMS\n", ""),
                "duplicate" => mutate(valid, "CMSCANVAS\t0.05\t1.1\n", "CMSCANVAS\t0.05\t1.1\nCMSCANVAS\t0.05\t1.1\n"),
                "after CANVAS" => mutate(mutate(valid, "CMSCANVAS\t0.05\t1.1\n", ""), r"(\nCANVAS\t[^\n]*\n)", s"\1CMSCANVAS\t0.05\t1.1\n"),
                "at the end" => mutate(valid, "CMSCANVAS\t0.05\t1.1\n", "") * "CMSCANVAS\t0.05\t1.1\n",
                "missing field" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\t0.05"),
                "extra field" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\t0.05\t1.1\t0"),
                "nan extra space" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\tnan\t1.1"),
                "infinite offset" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\t0.05\tinf"),
                "extra space range" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\t0.5\t1.1"),
                "negative extra space" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\t-0.01\t1.1"),
                "offset range" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\t0.05\t10.5"),
                "zero offset" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\t0.05\t0"),
                "text offset" => mutate(valid, "CMSCANVAS\t0.05\t1.1", "CMSCANVAS\t0.05\tdefault"),
                # Found by scripts/fuzz_native_requests.jl: ROOT's TAxis::SetTimeOffset crashed
                # (SIGSEGV in strftime) on a huge offset; offsets need |t| < 2^31.
                "huge time offset" => mutate(valid, request_line(valid, "AXIS\ty"),
                                             set_fields(request_line(valid, "AXIS\ty"), 24 => "1e308")),
                "time offset beyond 2^31" => mutate(valid, request_line(valid, "AXIS\tx"),
                                                    set_fields(request_line(valid, "AXIS\tx"), 24 => "2147483648")),
            ]
            for (name, broken) in cases
                @test broken != valid
                write(input, broken)
                ok = rejects_request(exe, input, joinpath(dir, "bad"))
                ok || @info "CMSCANVAS case not rejected" name
                @test ok
            end
            # A rejected stray row names the rule.
            write(input, cases[1].second)
            @test occursin("CMSCANVAS", run_native(exe, input, joinpath(dir, "bad")).stderr)
        end
    else
        @test_skip false
    end
end
