using ROOTPlots: stats_option, fit_option, stats_request, stats_box_position, write_hist_stats

"""Parse a ROOTPLOTS_STATS_DUMP file into boxes with printed lines and ROOT values."""
function read_stats_dump(path)
    boxes = NamedTuple[]
    current = nothing
    for line in readlines(path)
        f = split(line, '\t')
        if f[1] == "BOX"
            current = (geometry=parse.(Float64, f[2:5]), optstat=parse(Int, f[6]),
                       optfit=parse(Int, f[7]), lines=String[], values=Dict{String,Float64}(),
                       params=Tuple{String,Float64,Float64}[])
        elseif f[1] == "LINE"
            push!(current.lines, join(f[2:end], '\t'))
        elseif f[1] == "VALUE"
            current.values[f[2]] = parse(Float64, f[3])
        elseif f[1] == "PARAM"
            push!(current.params, (f[2], parse(Float64, f[3]), parse(Float64, f[4])))
        elseif f[1] == "END"
            push!(boxes, current)
        end
    end
    return boxes
end

"""Value printed after '=' on the box line whose label matches."""
function box_value(box, label)
    for line in box.lines
        m = match(r"^\s*(.*?)\s*=\s*(.*?)\s*$", line)
        m === nothing && continue
        m.captures[1] == label && return m.captures[2]
    end
    return nothing
end

function with_stats_dump(f)
    mktempdir() do dir
        path = joinpath(dir, "stats.txt")
        withenv("ROOTPLOTS_STATS_DUMP" => path) do
            f(dir)
        end
        return isfile(path) ? read_stats_dump(path) : NamedTuple[]
    end
end

"""FRAME records of a ROOTPLOTS_STATS_DUMP run: 1-D frame y range (pad coordinates),
highest data point under the boxes and lowest box edge (NDC)."""
function with_frame_dump(f)
    mktempdir() do dir
        path = joinpath(dir, "stats.txt")
        withenv("ROOTPLOTS_STATS_DUMP" => path) do
            f(dir)
        end
        frames = NamedTuple[]
        isfile(path) || return frames
        for line in readlines(path)
            v = split(line, '\t')
            v[1] == "FRAME" || continue
            push!(frames, (headroom=v[2] == "1", ymin=parse(Float64, v[3]), ymax=parse(Float64, v[4]),
                           top=parse(Float64, v[5]), box=parse(Float64, v[6])))
        end
        return frames
    end
end

@testset "statistics and fit box options" begin
    @test stats_option(true) == 1111
    @test stats_option(false) == 0 && stats_option(nothing) == 0
    @test stats_option((:name, :entries)) == 11
    @test stats_option([:entries, :mean_error, :std_error]) == 2210
    @test stats_option(:underflow) == 10000
    @test stats_option((:mean, :mean_error)) == 200      # the error variant wins
    @test stats_option((:integral, :overflow, :underflow)) == 1110000
    @test stats_option((:skewness_error, :kurtosis)) == 120000000
    @test stats_option("ksiourmen") == 111111111
    @test stats_option("KSIourRMen") == 222112211
    @test stats_option("") == 0
    @test stats_option(1) == 1111 && stats_option(1001111) == 1001111
    @test_throws ArgumentError stats_option(3)
    @test_throws ArgumentError stats_option(-1)
    @test_throws ArgumentError stats_option(1_000_000_000)
    @test_throws ArgumentError stats_option("x")
    @test_throws ArgumentError stats_option((:median,))
    @test_throws ArgumentError stats_option(1.0)
    @test_throws ArgumentError stats_option(())
    @test fit_option(true) == 1111 && fit_option(false) == 0
    @test fit_option((:chi2, :values)) == 101
    @test fit_option((:prob, :chi2, :errors, :values)) == 1111
    @test fit_option(:free_values) == 2
    @test fit_option("pcev") == 1111 && fit_option("ev") == 11
    @test fit_option(1) == 111 && fit_option(1011) == 1011
    @test_throws ArgumentError fit_option(10_000)
    @test_throws ArgumentError fit_option("m")
    @test_throws ArgumentError fit_option((:mean,))
    @test stats_box_position(:auto) == ("auto", (0., 0., 0., 0.))
    @test stats_box_position(:bottomleft)[1] == "bottomleft"
    @test stats_box_position((0.1, 0.2, 0.4, 0.5)) == ("manual", (0.1, 0.2, 0.4, 0.5))
    @test_throws ArgumentError stats_box_position(:middle)
    @test_throws ArgumentError stats_box_position((0.5, 0.2, 0.4, 0.5))
    @test_throws ArgumentError stats_box_position((0.1, 0.2, 1.4, 0.5))
    @test_throws ArgumentError stats_box_position((0.1, 0.2, 0.4))
    @test_throws ArgumentError stats_request(; fitstats=true)                 # no attached fit
    @test_throws ArgumentError stats_request(; stats=true, allow_stats=false)  # graphs/stacks
    @test !stats_request(; stats=true).headroom && stats_request(; stats=true, stats_headroom=true).headroom
    @test_throws ArgumentError stats_request(; stats=true, stats_headroom=1)
    @test_throws ArgumentError stats_request(; stats=true, stats_headroom=true, one_dimensional=false)
    @test_throws ArgumentError stats_request(; stats=true, stats_headroom=true, ylimits=(0.0, 1.0))
    r = stats_request(; stats=(:entries,), stats_position=:topleft, name="")
    @test r.optstat == 10 && r.optfit == 0 && r.corner == "topleft" && r.name == "Histogram"
end

@testset "histogram entries, flows, and unbinned moments" begin
    x = [-2.0, 0.25, 0.5, 1.5, 1.75, 3.0, 9.0]
    w = [1.0, 2.0, 1.0, 0.5, 1.0, 1.0, 3.0]
    r = histogram(x, [0.0, 1.0, 2.0, 3.0]; weights=w, threads=:off)
    h = r.hist
    @test h.entries == length(x)
    @test h.underflow == r.underflow == (1.0, 1.0)
    @test h.overflow == r.overflow == (4.0, 10.0)   # 3.0 is the upper edge: overflow
    inside = 0 .<= x .< 3
    m = h.moments
    @test m[1] ≈ sum(w[inside]) && m[2] ≈ sum(w[inside] .^ 2)
    @test m[3] ≈ sum(w[inside] .* x[inside]) && m[4] ≈ sum(w[inside] .* x[inside] .^ 2)
    # Chunked/threaded filling reduces moments in a fixed order.
    big = collect(range(-1.0, 4.0; length=300_000))
    a = histogram(big, [0.0, 1.0, 2.0, 3.0]; threads=:off).hist
    b = histogram(big, [0.0, 1.0, 2.0, 3.0]; threads=:on).hist
    @test a.moments == b.moments && a.entries == b.entries == length(big)
    # Directly constructed histograms keep backward-compatible defaults.
    plain = Hist1D([0.0, 1.0, 2.0], [1.0, 2.0], [1.0, 2.0])
    @test plain.entries === nothing && plain.moments === nothing
    @test plain.underflow == (0.0, 0.0) && plain.overflow == (0.0, 0.0)
    explicit = Hist1D([0.0, 1.0], [2.0], [2.0]; entries=5, underflow=(1, 1), overflow=(2, 2),
                      moments=(2, 2, 1, 0.75))
    @test explicit.entries == 5 && explicit.underflow == (1.0, 1.0) && explicit.moments[3] == 1.0
    @test_throws ArgumentError Hist1D([0.0, 1.0], [2.0], [2.0]; entries=-1)
    @test_throws ArgumentError Hist1D([0.0, 1.0], [2.0], [2.0]; entries=2.5)
    @test_throws ArgumentError Hist1D([0.0, 1.0], [2.0], [2.0]; underflow=(1.0, -1.0))
    @test_throws ArgumentError Hist1D([0.0, 1.0], [2.0], [2.0]; overflow=(NaN, 1.0))
    @test_throws ArgumentError Hist1D([0.0, 1.0], [2.0], [2.0]; moments=(1, 1, 1))
    @test_throws ArgumentError Hist1D([0.0, 1.0], [2.0], [2.0]; moments=(1, -1, 1, 1))
    # ROOT's Rebin keeps statistics; projections are plain visible-bin histograms.
    rebinned = rebin(h, [0.0, 2.0, 3.0])
    @test rebinned.entries == h.entries && rebinned.moments == h.moments
    @test rebinned.underflow == h.underflow && rebinned.overflow == h.overflow
    io = IOBuffer(); write_hist_stats(io, h)
    rows = split(chomp(String(take!(io))), '\n')
    @test rows[1] == join(("HSTATS", 7, 2, 4, h.moments...), '\t')
    @test rows[2] == "FLOW\t0\t1.0\t1.0" && rows[3] == "FLOW\t4\t4.0\t10.0"
    io = IOBuffer(); write_hist_stats(io, plain)
    @test startswith(String(take!(io)), "HSTATS\t-\t2\t0\n")

    xs, ys = [0.5, 1.5, -1.0, 0.5, 2.5], [0.5, 0.5, 0.5, 3.0, 1.5]
    r2 = histogram((xs, ys), ([0.0, 1.0, 2.0], [0.0, 1.0, 2.0]); threads=:off)
    h2 = r2.hist
    @test h2.entries == 5 && h2.bins_with_flow.sumw == r2.bins_with_flow.sumw
    @test h2.bins_with_flow.sumw !== r2.bins_with_flow.sumw # owned copy
    @test h2.moments == (2.0, 2.0, 2.0, 2.5, 1.0, 0.5, 1.0)
    io = IOBuffer(); write_hist_stats(io, h2)
    rows = split(chomp(String(take!(io))), '\n')
    @test rows[1] == join(("HSTATS", 5, 3, 7, h2.moments...), '\t')
    @test Set(rows[2:end]) == Set(["FLOW\t0\t1\t1.0\t1.0", "FLOW\t1\t3\t1.0\t1.0", "FLOW\t3\t2\t1.0\t1.0"])
    r3 = histogram(([0.5, 5.0], [0.5, 0.5], [0.5, 0.5]), ([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]))
    @test r3.hist.entries == 2 && length(r3.hist.moments) == 11 && r3.hist.moments[1] == 1.0
    @test projection(h2, :x).entries === nothing
    @test Hist2D([0., 1.], [0., 1.], fill(1., 1, 1), fill(1., 1, 1)).bins_with_flow === nothing
    @test_throws ArgumentError Hist2D([0., 1.], [0., 1.], fill(1., 1, 1), fill(1., 1, 1);
                                      bins_with_flow=(sumw=zeros(3, 3), sumw2=zeros(3, 3)))
    @test_throws ArgumentError Hist2D([0., 1.], [0., 1.], fill(1., 1, 1), fill(1., 1, 1); moments=(1, 1))
end

@testset "statistics box keyword validation (no ROOT needed)" begin
    mktempdir() do dir
        h = histogram([0.5, 1.5], [0.0, 1.0, 2.0]).hist
        spec = GraphSpec(legend=:none)
        g = GraphSeries([1.0, 2.0, 3.0], [1.0, 2.0, 3.1]; yerrors=[0.1, 0.1, 0.1])
        bad = FitResult(:pol0, (0.0, 2.0), [1.0], [0.1], fill(0.01, 1, 1), falses(1), 4, false,
                        0, 1.0, :chi2, 1, 2, 0.0, "Minuit2", false)
        stem = joinpath(dir, "never")
        @test_throws ArgumentError render_root(h, spec, stem; draw="E1", fitstats=true)
        @test_throws ArgumentError render_root(h, spec, stem; draw="E1", stats=:median)
        @test_throws ArgumentError render_root(h, spec, stem; draw="E1", stats=true, stats_position=:center)
        @test_throws ArgumentError render_root(h, spec, stem; draw=:hist, attached_fit=bad)
        @test_throws ArgumentError render_root([g], spec, stem; fitstats=true)
        @test_throws ArgumentError render_root([g], spec, stem; attached_fit=bad)
        @test_throws ArgumentError render_root([g], spec, stem; attached_to=2)
        @test_throws ArgumentError render_root(HistStack([h, h]), spec, stem; stats=true)
        @test_throws ArgumentError render_root(PlotPanel([g], spec; stats=true), stem)
        @test_throws ArgumentError render_root(PlotPanel(h, spec; stat=true), stem)
        @test_throws ArgumentError render_root([h, h], spec, stem; draw=["HIST", "HIST SAMES"], stats="q")
        # stats_headroom: Bool only, 1-D frames only, not with explicit ylimits.
        fixed = GraphSpec(legend=:none, ylimits=(0.0, 5.0))
        h2 = histogram(([0.5], [0.5]), ([0.0, 1.0], [0.0, 1.0])).hist
        profile = Profile1D([0.0, 1.0, 2.0], [0.5, 1.5], [1.0, 2.0])
        @test_throws ArgumentError render_root(h, spec, stem; draw="E1", stats=true, stats_headroom=1)
        @test_throws ArgumentError render_root(h, fixed, stem; draw="E1", stats=true, stats_headroom=true)
        @test_throws ArgumentError render_root(h, fixed, stem; draw=:hist, stats=true, stats_headroom=true)
        @test_throws ArgumentError render_root([g], fixed, stem; attached_fit=nothing, stats_headroom=true)
        @test_throws ArgumentError render_root(profile, fixed, stem; stats=true, stats_headroom=true)
        @test_throws ArgumentError render_root(h2, Graph2DSpec(), stem; stats=true, stats_headroom=true)
        @test_throws ArgumentError render_root(PlotPanel(h, fixed; stats=true, stats_headroom=true), stem)
        @test !isfile(stem * ".png")
    end
end

if TEST_NATIVE_ROOT
    @testset "native statistics and fit boxes" begin
        # Deterministic, roughly Gaussian sample (Irwin-Hall of four low-discrepancy sequences).
        u(i, a) = mod(i * a, 1.0)
        x = vcat([(u(i, 0.6180339887) + u(i, 0.4142135623) + u(i, 0.7320508075) +
                   u(i, 0.2360679774) - 2) * sqrt(3.0) for i in 1:2000], fill(-5.0, 3), fill(6.0, 2))
        edges = collect(range(-3.0, 3.0; length=25))
        r = histogram(x, edges)
        h = r.hist
        inside = filter(v -> -3 <= v < 3, x)
        exact_mean = sum(inside) / length(inside)
        exact_std = sqrt(sum(inside .^ 2) / length(inside) - exact_mean^2)
        spec = GraphSpec(xlabel="x", ylabel="Entries")
        fit = fit_root(h, :gaus; initial=[100.0, 0.0, 1.0])
        @test fit.status == 0 && fit.valid
        boxes = with_stats_dump() do dir
            render_root(h, spec, joinpath(dir, "default"); draw="E1", stats=true, label="Sample")
            render_root(h, spec, joinpath(dir, "steps"); stats=true)   # :steps becomes HIST
            render_root(h, spec, joinpath(dir, "fields"); draw=:hist,
                        stats=(:entries, :underflow, :overflow, :integral), stats_position=:topleft)
            render_root(h, spec, joinpath(dir, "fit"); draw=:e1, stats=true, fitstats=true,
                        attached_fit=fit, fit_label="Gaussian")
            render_root(h, spec, joinpath(dir, "fitonly"); draw="E1", fitstats=(:chi2, :values),
                        attached_fit=fit, stats_position=(0.15, 0.6, 0.45, 0.85))
            render_root(Hist1D(h.edges, h.sumw, h.sumw2), spec, joinpath(dir, "plain"); draw="E1", stats=true)
            render_root(h, spec, joinpath(dir, "nobox"); draw="E1")   # no dump entry
        end
        @test length(boxes) == 6
        b = boxes[1]
        @test b.optstat == 1111 && b.optfit == 0
        @test b.lines[1] == "Sample"
        @test parse(Int, box_value(b, "Entries")) == length(x) == b.values["entries"]
        @test b.values["mean1"] ≈ exact_mean rtol=1e-12
        @test b.values["std1"] ≈ exact_std rtol=1e-12
        @test parse(Float64, box_value(b, "Mean")) ≈ exact_mean rtol=1e-3
        @test parse(Float64, box_value(b, "Std Dev")) ≈ exact_std rtol=1e-3
        @test b.geometry ≈ [0.78, 0.775, 0.98, 0.935] atol=1e-6   # ROOT's default placement
        @test boxes[2].lines[1] == "Histogram" && boxes[2].values["entries"] == length(x)
        f = boxes[3]
        @test f.optstat == 1110010 && f.geometry[1] < 0.2
        @test parse(Int, box_value(f, "Underflow")) == count(<(-3), x)
        @test parse(Int, box_value(f, "Overflow")) == count(>=(3), x)
        @test parse(Float64, box_value(f, "Integral")) == length(inside)
        fb = boxes[4]
        @test fb.optstat == 1111 && fb.optfit == 1111
        @test [p[1] for p in fb.params] == ["Constant", "Mean", "Sigma"]
        @test [p[2] for p in fb.params] ≈ fit.parameters rtol=1e-12
        @test [p[3] for p in fb.params] ≈ fit.errors rtol=1e-12
        @test fb.values["chi2"] ≈ fit.chi2 rtol=1e-12
        @test fb.values["ndf"] == fit.ndf
        @test fb.values["prob"] ≈ chi2_pvalue(fit) rtol=1e-8
        @test any(l -> startswith(l, "#chi^{2} / ndf"), fb.lines) && box_value(fb, "Prob") !== nothing
        fo = boxes[5]
        @test fo.optstat == 0 && fo.optfit == 101 && length(fo.lines) == 4
        @test fo.geometry ≈ [0.15, 0.6, 0.45, 0.85] atol=1e-6
        p = boxes[6]   # entries unknown: ROOT's effective entries and bin-centre moments
        @test p.values["entries"] ≈ effective_entries(Hist1D(h.edges, h.sumw, h.sumw2))
        @test p.values["mean1"] ≈ binned_mean(h) rtol=1e-12

        h2 = histogram(([0.5, 1.5, 1.5, 2.5, -1.0], [0.5, 0.5, 1.5, 1.5, 0.5]),
                       ([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0])).hist
        profile = Profile1D([0.0, 1.0, 2.0], [0.5, 0.5, 1.5, 2.5], [1.0, 3.0, 5.0, 7.0])
        g = GraphSeries([1.0, 2.0, 3.0, 4.0], [2.1, 3.9, 6.2, 7.8]; yerrors=fill(0.2, 4))
        gfit = fit_root(g, :pol1; initial=[0.0, 2.0])
        nd = with_stats_dump() do dir
            render_root(h2, Graph2DSpec(), joinpath(dir, "h2"); stats=true, label="map")
            render_root(profile, GraphSpec(), joinpath(dir, "profile"); stats=true)
            render_root([g], GraphSpec(legend=:none), joinpath(dir, "graph"); attached_fit=gfit,
                        fitstats=true, stats_position=:bottomright)
            render_root([h, histogram(x .+ 0.5, edges).hist], spec, joinpath(dir, "overlay");
                        draw=["HIST", "HIST SAMES"], labels=["A", "B"], stats=true)
            render_root([h, h], spec, joinpath(dir, "overlay_same"); draw=["HIST", "HIST SAME"], stats=true)
            render_root(PlotFigure([PlotPanel(h, spec; draw=:e1, stats=true, attached_fit=fit, fitstats=true),
                                    PlotPanel(h2, Graph2DSpec(); stats=(:entries, :mean)),
                                    PlotPanel(profile, GraphSpec(); stats=true),
                                    PlotPanel([g], GraphSpec(); attached_fit=gfit, fitstats=true)];
                                   layout=(2, 2), size=(800, 600)), joinpath(dir, "figure"))
        end
        @test length(nd) == 2 + 1 + 2 + 1 + 4
        m = nd[1]
        @test m.lines[1] == "map" && m.values["entries"] == 5
        @test m.values["mean1"] ≈ 1.5 && m.values["mean2"] ≈ 1.0
        @test box_value(m, "Mean x") !== nothing && box_value(m, "Std Dev y") !== nothing
        @test m.geometry[3] < 1 - 0.18 + 1e-6    # moved left of the color palette
        pr = nd[2]
        @test pr.lines[1] == "Profile" && pr.values["entries"] == 4 && box_value(pr, "Mean y") !== nothing
        gb = nd[3]
        @test gb.optstat == 0 && [q[1] for q in gb.params] == ["p0", "p1"]
        @test [q[2] for q in gb.params] ≈ gfit.parameters rtol=1e-12
        @test gb.geometry[2] < 0.2               # bottom-right corner
        a, bb = nd[4], nd[5]
        @test a.lines[1] == "A" && bb.lines[1] == "B"
        @test bb.geometry[4] ≈ a.geometry[2] atol=1e-6   # stacked directly below
        @test nd[6].lines[1] == "Histogram"              # SAME layer: no second box
        fig = nd[7:10]
        @test fig[1].optfit == 1111 && length(fig[1].params) == 3
        @test fig[2].optstat == 110 && length(fig[2].lines) == 3
        @test fig[3].lines[1] == "Profile"
        @test fig[4].optstat == 0 && length(fig[4].params) == 2
    end

    @testset "statistics box headroom" begin
        u(i, a) = mod(i * a, 1.0)
        x = [(u(i, 0.6180339887) + u(i, 0.4142135623) + u(i, 0.7320508075) +
              u(i, 0.2360679774) - 2) * sqrt(3.0) + 0.5 for i in 1:3000]
        h = histogram(x, collect(range(-3.0, 3.0; length=31))).hist
        fit = fit_root(h, :gaus; initial=[150.0, 0.5, 1.0])
        @test fit.status == 0 && fit.valid
        spec = GraphSpec(xlabel="x", ylabel="Entries", legend=:topleft)
        px = [4 * u(i, 0.6180339887) for i in 1:500]
        profile = Profile1D(collect(0.0:0.5:4.0), px, 2 .* px .+ 0.1 .* sin.(7 .* px))
        g = GraphSeries(collect(1.0:6.0), [1.1, 2.0, 2.9, 4.2, 5.0, 6.1]; yerrors=fill(0.2, 6))
        gfit = fit_root(g, :pol1; initial=[0.0, 1.0])
        box = (stats=true, fitstats=true, attached_fit=fit)
        frames = with_frame_dump() do dir
            render_root(h, spec, joinpath(dir, "plain"); draw=:e1, box...)                        # 1
            render_root(h, spec, joinpath(dir, "headroom"); draw=:e1, box..., stats_headroom=true) # 2
            render_root(h, spec, joinpath(dir, "raw"); draw="HIST", stats=true, stats_headroom=true) # 3
            render_root(PlotFigure([PlotPanel(h, spec; draw=:e1, box..., stats_headroom=true),     # 4
                                    PlotPanel(profile, GraphSpec(legend=:none); stats=true,      # 5
                                              stats_headroom=true)]; layout=(1, 2), size=(900, 400)),
                        joinpath(dir, "figure"))
            render_root(profile, GraphSpec(legend=:none), joinpath(dir, "profile"); stats=true,
                        stats_headroom=true)                                                     # 6
            render_root([g], GraphSpec(legend=:none), joinpath(dir, "graph"); attached_fit=gfit,
                        fitstats=true, stats_headroom=true)                                      # 7
            for headroom in (false, true)                                                        # 8, 9
                render_root(profile, GraphSpec(legend=:none), joinpath(dir, "low$headroom"); stats=true,
                            stats_position=:bottomright, stats_headroom=headroom)
            end
        end
        @test length(frames) == 9
        plain, raised = frames[1], frames[2]
        @test !plain.headroom && plain.top > plain.box          # ROOT's default box covers data
        @test raised.headroom && raised.top <= raised.box - 0.02 + 1e-6
        @test raised.ymin ≈ plain.ymin && raised.ymax > plain.ymax
        for f in frames[3:7]
            @test f.headroom && f.top <= f.box - 0.02 + 1e-6
        end
        # A bottom box cannot be cleared by a higher maximum: the frame is unchanged.
        @test frames[8].ymin == frames[9].ymin && frames[8].ymax == frames[9].ymax
    end

    @testset "profile error bars use one colour" begin
        u(i, a) = mod(i * a, 1.0)
        px = [4 * u(i, 0.6180339887) for i in 1:2000]
        py = 2 .* px .+ 0.5 .* [u(i, 0.4142135623) - 0.5 for i in 1:2000]
        profile = Profile1D(collect(0.0:0.25:4.0), px, py)
        mktempdir() do dir
            # ROOT's default TProfile attributes painted every other E1 bar blue.
            single = render_root(profile, GraphSpec(legend=:none), joinpath(dir, "profile")).png
            @test ROOTPlots.image_colored_fraction(single) == 0
            figure = render_root(PlotFigure([PlotPanel(profile, GraphSpec(legend=:none)),
                                             PlotPanel(profile, GraphSpec(legend=:none); stats=true)];
                                            layout=(1, 2), size=(800, 400)), joinpath(dir, "figure")).png
            @test ROOTPlots.image_colored_fraction(figure) == 0
            # Positive control: the measure sees coloured error bars.
            h = Hist1D([0.0, 1.0, 2.0], [3.0, 4.0], [3.0, 4.0])
            red = render_root(h, GraphSpec(legend=:none), joinpath(dir, "red"); draw=:e1,
                              style=GraphStyle(color="#cc0000")).png
            @test ROOTPlots.image_colored_fraction(red) > 0
        end
    end
end
