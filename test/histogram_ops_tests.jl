# Scalar scaling, merging, and flow-aware projections (src/histogram_ops.jl).
# Every operation is checked against a direct fill of the equivalent events;
# expected numbers marked "ROOT 6.40.04" were printed by a plain ROOT macro
# (TH1::Scale, TH1::Add/Merge, TH2::ProjectionX/Y, TH3::Project3D).

"""All numbers of two histograms agree (rtol for floating sums); entries exactly."""
function same_hist(a::Hist1D, b::Hist1D; rtol=1e-12)
    a.edges == b.edges && a.entries == b.entries &&
        isapprox([a.underflow[1], a.sumw..., a.overflow[1]], [b.underflow[1], b.sumw..., b.overflow[1]]; rtol, atol=1e-12) &&
        isapprox([a.underflow[2], a.sumw2..., a.overflow[2]], [b.underflow[2], b.sumw2..., b.overflow[2]]; rtol, atol=1e-12) &&
        same_moments(a.moments, b.moments, rtol)
end

function same_hist(a::T, b::T; rtol=1e-12) where {T<:Union{Hist2D,Hist3D}}
    ROOTPlots.hist_nd_edges(a) == ROOTPlots.hist_nd_edges(b) && a.entries == b.entries &&
        isapprox(a.sumw, b.sumw; rtol, atol=1e-12) && isapprox(a.sumw2, b.sumw2; rtol, atol=1e-12) &&
        (a.bins_with_flow === nothing) == (b.bins_with_flow === nothing) &&
        (a.bins_with_flow === nothing ||
         (isapprox(a.bins_with_flow.sumw, b.bins_with_flow.sumw; rtol, atol=1e-12) &&
          isapprox(a.bins_with_flow.sumw2, b.bins_with_flow.sumw2; rtol, atol=1e-12))) &&
        same_moments(a.moments, b.moments, rtol)
end

same_moments(a, b, rtol) = (a === nothing && b === nothing) ||
    (a !== nothing && b !== nothing && isapprox(collect(a), collect(b); rtol, atol=1e-12))

"""The visible part of an ND result equals its stored full array exactly."""
nd_consistent(h) = h.bins_with_flow === nothing ||
    (v = map(e -> 2:length(e), ROOTPlots.hist_nd_edges(h));
     h.bins_with_flow.sumw[v...] == h.sumw && h.bins_with_flow.sumw2[v...] == h.sumw2)

"""No array of one histogram is an array of the other (results never alias inputs)."""
no_alias(a, b) = !any(x -> any(y -> x === y, ROOTPlots.data_arrays(b)), ROOTPlots.data_arrays(a)) &&
                 !any(x -> any(y -> x === y, ROOTPlots.hist_edges(b)), ROOTPlots.hist_edges(a))

in_axis(c, e) = e[1] <= c < e[end]

"""True when no nonzero-weight event is visible on the kept axes but in a flow bin of another."""
projection_clean(coords, edges, w, kept) = !any(eachindex(w)) do k
    w[k] != 0 && all(j -> in_axis(coords[j][k], edges[j]), kept) &&
        any(j -> !(j in kept) && !in_axis(coords[j][k], edges[j]), eachindex(coords))
end

# Events with in-range, flow, signed, and zero weights.
const OPS_X = [-1.0, 0.5, 1.5, 1.5, 2.5, 3.5, 5.0, 0.2, 3.9]
const OPS_W = [2.0, 1.0, -0.5, 3.0, 2.0, 1.0, 4.0, 0.0, -1.25]
const OPS_Y = [0.5, 2.5, 1.5, -2.0, 2.5, 0.5, 1.5, 2.9, 7.0]
const OPS_Z = [0.5, 1.5, 0.5, 1.5, 3.0, 0.5, -0.1, 1.0, 1.9]
const OPS_E = ([0.0, 1, 2, 3, 4], [0.0, 1, 2, 3], [0.0, 1, 2])
const OPS_COLS = ((OPS_X,), (OPS_X, OPS_Y), (OPS_X, OPS_Y, OPS_Z))

ops_fill(cols::Tuple{Any}, w) = histogram(cols[1], OPS_E[1]; weights=w, threads=:off).hist
ops_fill(cols::Tuple, w) = histogram(cols, OPS_E[1:length(cols)]; weights=w, threads=:off).hist

"""A copy of the histogram with another entry count."""
with_entries(h::Hist1D, n) = Hist1D(h.edges, h.sumw, h.sumw2; underflow=h.underflow, overflow=h.overflow,
                                    entries=n, moments=h.moments)
with_entries(h, n) = typeof(h)(ROOTPlots.hist_nd_edges(h)..., h.sumw, h.sumw2;
                               bins_with_flow=h.bins_with_flow, entries=n, moments=h.moments)

"""A copy of a Hist1D without statistics (entries and moments unknown)."""
no_stats(h::Hist1D) = Hist1D(h.edges, h.sumw, h.sumw2; underflow=h.underflow, overflow=h.overflow)

@testset "histogram scaling equals a direct fill with scaled weights" begin
    for c in (2.5, -2.5, 1/3, 0, -0.0, 3//4, 1e-3), cols in OPS_COLS
        direct = ops_fill(cols, Float64(c) .* OPS_W)
        scaled = isodd(length(cols)) ? c * ops_fill(cols, OPS_W) : ops_fill(cols, OPS_W) * c
        @test same_hist(direct, scaled)
        @test scaled.entries == length(OPS_X)
        @test scaled isa Hist1D || nd_consistent(scaled)
    end
    h = ops_fill((OPS_X,), OPS_W)
    # ROOT 6.40.04, the first seven events after TH1::Scale(-2.5).
    r = -2.5 * ops_fill((OPS_X[1:7],), OPS_W[1:7])
    @test r.sumw == [-2.5, -6.25, -5.0, -2.5] && r.sumw2 == [6.25, 57.8125, 25.0, 6.25]
    @test r.underflow == (-5.0, 25.0) && r.overflow == (-10.0, 100.0)
    @test r.moments == (-16.25, 95.3125, -31.875, -76.5625) && r.entries == 7
    # Scale(0): zero contents and statistics, fill count kept, no negative zeros.
    z = 0 * h
    @test all(iszero, z.sumw) && all(iszero, z.sumw2) && z.underflow == (0.0, 0.0)
    @test z.moments == (0.0, 0.0, 0.0, 0.0) && z.entries == h.entries
    @test !any(signbit, z.sumw) && !any(signbit, z.moments) && !any(signbit, (-0.0 * h).sumw)
    # h / c is h * inv(c), like ROOT Scale(1/c).
    @test same_hist(h / 4, 0.25 * h; rtol=0) && same_hist(h / -0.5, -2 * h; rtol=0)
    # Unknown statistics stay unknown; ND histograms without flow arrays keep none.
    u = Hist1D([0, 1, 2], [1, -2], [1, 4])
    @test (3u).entries === nothing && (3u).moments === nothing && (3u).sumw == [3, -6] && (3u).sumw2 == [9, 36]
    u2 = Hist2D([0, 1, 2], [0, 1], reshape([1.0, 2], 2, 1), reshape([1.0, 4], 2, 1); entries=5)
    @test (2u2).bins_with_flow === nothing && (2u2).entries == 5 && vec((2u2).sumw2) == [4, 16]
    # No mutation and no aliasing.
    h3 = ops_fill(OPS_COLS[3], OPS_W)
    before = deepcopy(h3)
    s3 = 2 * h3
    @test same_hist(h3, before; rtol=0) && no_alias(s3, h3)
    s3.bins_with_flow.sumw[2, 2, 2] = 99.0
    s3.sumw[1, 1, 1] = 99.0
    @test same_hist(h3, before; rtol=0)
    @test no_alias(2h, h) && no_alias(h * 1, h) && no_alias(h / 1, h)
    # Errors: non-finite factor, division by zero, Float64 overflow (contents or c²).
    for bad in (Inf, -Inf, NaN)
        @test_throws ArgumentError bad * h
        @test_throws ArgumentError h / bad
    end
    @test_throws ArgumentError h / 0
    @test_throws ArgumentError h / 0.0
    @test_throws ArgumentError 1e200 * h
    huge = Hist1D([0, 1], [1e300], [1.0])
    @test_throws ArgumentError 1e10 * huge
    @test_throws ArgumentError 1e300 * Hist1D([0, 1], [1.0], [1.0])
    @test (1e-300 * huge).sumw == [1.0]
    # Arrays modified after construction are revalidated before use.
    damaged = ops_fill(OPS_COLS[2], OPS_W)
    damaged.sumw[1, 1] += 1
    @test_throws ArgumentError 2 * damaged
    @test_throws ArgumentError merge(damaged, ops_fill(OPS_COLS[2], OPS_W))
    @test_throws ArgumentError projection(damaged, :x; flow=true)
end

@testset "merge equals a direct fill of independent samples" begin
    a, b = 1:5, 6:9
    for cols in OPS_COLS
        direct = ops_fill(cols, OPS_W)
        ha, hb = ops_fill(map(c -> c[a], cols), OPS_W[a]), ops_fill(map(c -> c[b], cols), OPS_W[b])
        before = deepcopy(ha)
        m = merge(ha, hb)
        @test same_hist(direct, m) && m.entries == length(OPS_X)
        @test same_hist(ha, before; rtol=0) && no_alias(m, ha) && no_alias(m, hb)
        @test m isa Hist1D || nd_consistent(m)
        # Three inputs, one with zero weights only: it adds entries, nothing else.
        hz = ops_fill(map(c -> c[1:2], cols), [0.0, 0.0])
        m3 = merge(hz, ha, hb)
        @test m3.entries == length(OPS_X) + 2 && same_hist(with_entries(m3, length(OPS_X)), direct)
        @test_throws ArgumentError merge(ha, ha)
        @test_throws ArgumentError merge(ha, hb, ha)
    end
    # ROOT 6.40.04: the first seven events plus three more via TH1::Add (and TH1::Merge).
    h = ops_fill((OPS_X[1:7],), OPS_W[1:7])
    c = ops_fill(([0.5, -3.0, 3.2],), [2.0, 1.0, -1.0])
    m = merge(h, c)
    @test m.entries == 10 && isapprox(collect(m.moments), [7.5, 20.25, 10.55, 20.885]; rtol=1e-12)
    # Unknown entries or moments are never mixed with known ones.
    known = Hist1D([0, 1, 2], [1, 2], [1, 2]; entries=3, moments=(3.0, 3.0, 4.0, 6.0))
    nomom = Hist1D([0, 1, 2], [1, 2], [1, 2]; entries=3)
    noent = Hist1D([0, 1, 2], [1, 2], [1, 2]; moments=(3.0, 3.0, 4.0, 6.0))
    @test merge(known, 2 * known).entries == 6 && merge(known, 2 * known).moments == (9.0, 15.0, 12.0, 18.0)
    @test merge(known, nomom).entries == 6 && merge(known, nomom).moments === nothing
    @test merge(known, noent).entries === nothing && merge(known, noent).moments == (6.0, 6.0, 8.0, 12.0)
    # ND: flow arrays from one input, none (zero flow bins) in the other.
    f2 = ops_fill(OPS_COLS[2], OPS_W)
    plain = Hist2D(f2.xedges, f2.yedges, f2.sumw, f2.sumw2)
    mp = merge(plain, f2)
    @test mp.bins_with_flow !== nothing && nd_consistent(mp) && mp.entries === nothing && mp.moments === nothing
    @test mp.sumw == 2 .* f2.sumw && mp.bins_with_flow.sumw[1, :] == f2.bins_with_flow.sumw[1, :]
    @test merge(plain, 2 * plain).bins_with_flow === nothing
    # Edges, types, overflow, and the same sample twice.
    @test_throws ArgumentError merge(Hist1D([0, 1, 2], [1, 2], [1, 2]), Hist1D([0, 1, 3], [1, 2], [1, 2]))
    @test_throws ArgumentError merge(f2, ops_fill((OPS_X,), OPS_W))
    @test_throws ArgumentError merge(f2, Hist2D([0, 1, 2, 3, 5], [0, 1, 2, 3], f2.sumw, f2.sumw2))
    err = try merge(h, h); nothing catch e; e end
    @test err isa ArgumentError && occursin("2 * h", err.msg)
    err = try merge(h, Hist1D([0, 2, 4], [1, 1], [1, 1])); nothing catch e; e end
    @test err isa ArgumentError && occursin("rebin", err.msg)
    @test_throws ArgumentError merge(Hist1D([0, 1], [1e308], [1.0]), Hist1D([0, 1], [1e308], [1.0]))
    @test_throws ArgumentError merge(Hist1D([0, 1], [1.0], [1.0]; entries=typemax(Int)),
                                     Hist1D([0, 1], [1.0], [1.0]; entries=1))
end

@testset "projection flow policy" begin
    e2 = ([0.0, 1, 2, 3], [0.0, 1, 2, 3])
    # Default: visible bins only, statistics unknown (unchanged behaviour).
    h2 = ops_fill(OPS_COLS[2], OPS_W)
    p = projection(h2, :x)
    @test p.sumw ≈ vec(sum(h2.sumw; dims=2)) && p.underflow == (0.0, 0.0) && p.entries === nothing && p.moments === nothing
    # flow=true: all events are in the result; here some hit integrated-axis flow bins.
    for (k, axis) in enumerate((:x, :y))
        pf = projection(h2, axis; flow=true)
        d = histogram(OPS_COLS[2][k], OPS_E[k]; weights=OPS_W, threads=:off).hist
        @test same_hist(pf, no_stats(d))
    end
    # Safe case: no event with visible x has a flow y; projected-axis (x) flows are allowed.
    xs = [0.5, 1.5, 2.5, 2.5, 5.0, -1.0, -2.0, 0.2]
    ys = [0.5, 2.5, 1.5, 2.5, 1.5, 9.0, -1.0, 0.0]
    ws = [1.0, 2.0, 0.5, 1.5, 3.0, -1.0, 2.0, 0.0]
    s2 = histogram((xs, ys), e2; weights=ws, threads=:off).hist
    px = projection(s2, :x; flow=true)
    @test same_hist(px, histogram(xs, e2[1]; weights=ws, threads=:off).hist) && px.entries == length(xs)
    @test no_alias(px, s2)
    # ROOT 6.40.04 ("inside": the first four events): ProjectionX/Y reuse entries 4 and
    # (Σw, Σw², Σwx, Σwx²) = (5, 7.5, 8.5, 17.25) / (5, 7.5, 10, 23.25).
    inside = histogram((xs[1:4], ys[1:4]), e2; weights=ws[1:4], threads=:off).hist
    @test projection(inside, :x; flow=true).moments == (5.0, 7.5, 8.5, 17.25)
    @test projection(inside, :y; flow=true).moments == (5.0, 7.5, 10.0, 23.25)
    @test projection(inside, :x; flow=true).entries == 4
    # ROOT 6.40.04 "projflow" (plus an x-overflow event): ROOT falls back to bin-centre
    # moments (equal here: the events sit at bin centres) and effective entries 3.333...;
    # here the exact moments and the fill count 5 are kept.
    pf = projection(histogram(([xs[1:4]; 5.0], [ys[1:4]; 1.5]), e2; weights=[ws[1:4]; 3.0], threads=:off).hist,
                    :x; flow=true)
    @test pf.moments == (5.0, 7.5, 8.5, 17.25) && pf.entries == 5 && pf.overflow == (3.0, 9.0)
    # Contaminated: one more event with visible x and overflowing y.
    c2 = histogram(([xs; 1.5], [ys; 5.0]), e2; weights=[ws; 3.0], threads=:off).hist
    pc = projection(c2, :x; flow=true)
    @test pc.entries === nothing && pc.moments === nothing && pc.sumw[2] == px.sumw[2] + 3
    # Signed cancellation in a flow cell (Σw = 0, Σw² = 8): conservative, statistics unknown.
    # (ROOT 6.40.04 "cancel" reuses them because the total content still matches.)
    k2 = histogram(([xs; 1.5; 1.5], [ys; 5.0; 6.0]), e2; weights=[ws; 2.0; -2.0], threads=:off).hist
    pk = projection(k2, :x; flow=true)
    @test pk.entries === nothing && pk.moments === nothing && pk.sumw2[2] == px.sumw2[2] + 8
    # A zero-weight event in an integrated-axis flow cell leaves the projection exact.
    z2 = histogram(([xs; 1.5], [ys; 5.0]), e2; weights=[ws; 0.0], threads=:off).hist
    pz = projection(z2, :x; flow=true)
    @test pz.entries == length(xs) + 1 && pz.moments == px.moments
    # Direct-constructor ND histograms (no flow arrays): flow=true sums the same visible
    # bins and keeps the given statistics.
    given = Hist2D(inside.xedges, inside.yedges, inside.sumw, inside.sumw2; entries=4, moments=inside.moments)
    pg = projection(given, :y; flow=true)
    @test pg.sumw == projection(given, :y).sumw && pg.entries == 4 && pg.moments == (5.0, 7.5, 10.0, 23.25)
    # 3D -> 1D and 3D -> 2D in every axis order, against direct fills of the projected coordinates.
    x3 = [0.5, 1.5, 2.5, 2.5, 0.2, 7.0, 1.1]
    y3 = [0.5, 2.5, 1.5, 2.5, 1.7, 1.5, -3.0]
    z3 = [0.5, 1.5, 0.5, 1.5, 0.9, 0.5, 9.0]
    w3 = [1.0, 2.0, 0.5, 1.5, -0.75, 3.0, 2.0]
    e3 = ([0.0, 1, 2, 3], [0.0, 1, 2, 3], [0.0, 1, 2])
    coords = (x3, y3, z3)
    h3 = histogram(coords, e3; weights=w3, threads=:off).hist
    names = (:x, :y, :z)
    for i in 1:3     # event 7 (x visible, y/z flow) and event 6 (x overflow) contaminate all three
        q = projection(h3, names[i]; flow=true)
        d = histogram(coords[i], e3[i]; weights=w3, threads=:off).hist
        @test !projection_clean(coords, e3, w3, (i,))
        @test same_hist(q, no_stats(d))
    end
    expected_clean = Dict((1, 2) => true, (2, 1) => true, (1, 3) => true, (3, 1) => true,
                          (2, 3) => false, (3, 2) => false)
    for ((i, j), clean) in expected_clean
        q = projection(h3, (names[i], names[j]); flow=true)
        d = histogram((coords[i], coords[j]), (e3[i], e3[j]); weights=w3, threads=:off).hist
        @test projection_clean(coords, e3, w3, (i, j)) == clean
        @test nd_consistent(q) && (q.xedges, q.yedges) == (e3[i], e3[j]) && no_alias(q, h3)
        if clean   # exact: equal to the direct 2D fill, entries and cross term included
            @test same_hist(q, d)
        else
            @test q.entries === nothing && q.moments === nothing
            @test same_hist(q, Hist2D(d.xedges, d.yedges, d.sumw, d.sumw2; bins_with_flow=d.bins_with_flow))
        end
    end
    # ROOT 6.40.04 ("3D 0", four in-range events): Project3D("xz") has z on X and x on Y with
    # stats (5, 7.5, 6, 8.25, 8.5, 17.25, 11); Project3D("yx") (5, 7.5, 8.5, 17.25, 10, 23.25, 19).
    t3 = histogram(map(c -> c[1:4], coords), e3; weights=w3[1:4], threads=:off).hist
    @test projection(t3, (:z, :x); flow=true).moments == (5.0, 7.5, 6.0, 8.25, 8.5, 17.25, 11.0)
    @test projection(t3, (:x, :y); flow=true).moments == (5.0, 7.5, 8.5, 17.25, 10.0, 23.25, 19.0)
    @test projection(t3, (:z, :x); flow=true).entries == 4
    @test projection(t3, (:z, :x)).entries === nothing    # default: unchanged, no statistics
    @test_throws ArgumentError projection(h3, (:x, :x); flow=true)
    @test_throws ArgumentError projection(h2, :z; flow=true)
end

@testset "ROOT statistics boxes after scaling, merging, and projection" begin
    if TEST_NATIVE_ROOT
        h = ops_fill((OPS_X[1:7],), OPS_W[1:7])
        c = ops_fill(([0.5, -3.0, 3.2],), [2.0, 1.0, -1.0])
        e2 = ([0.0, 1, 2, 3], [0.0, 1, 2, 3])
        xs, ys, ws = [0.5, 1.5, 2.5, 2.5], [0.5, 2.5, 1.5, 2.5], [1.0, 2.0, 0.5, 1.5]
        inside = histogram((xs, ys), e2; weights=ws, threads=:off).hist
        intflow = histogram(([xs; 1.5], [ys; 5.0]), e2; weights=[ws; 3.0], threads=:off).hist
        both = histogram(([xs; 1.5; 5.0], [ys; 5.0; 1.5]), e2; weights=[ws; 3.0; 4.0], threads=:off).hist
        cases = [-2.5 * h, merge(h, c), 0 * h, projection(inside, :x; flow=true),
                 projection(intflow, :x; flow=true), projection(both, :x; flow=true)]
        boxes = with_stats_dump() do dir
            for (k, x) in enumerate(cases)
                render_root(x, GraphSpec(), joinpath(dir, "op$k"); draw=:hist, stats=(:entries, :mean, :std))
            end
        end
        @test length(boxes) == length(cases)
        # Values printed by a plain ROOT 6.40.04 macro for the same fills:
        # Scale(-2.5), Add, Scale(0), ProjectionX of "inside" and of "intflow".
        expected = [(7, 1.961538462, 0.9294650749), (10, 1.406666667, 0.8977502746), (7, 0.0, 0.0),
                    (4, 1.7, sqrt(17.25 / 5 - 1.7^2)), (3.878787879, 1.625, 0.5994789404)]
        for (box, (n, mean, std)) in zip(boxes, expected)
            @test isapprox(box.values["entries"], n; rtol=1e-9)
            @test isapprox(box.values["mean1"], mean; rtol=1e-9, atol=1e-12)
            @test isapprox(box.values["std1"], std; rtol=1e-9, atol=1e-12)
        end
        # Fallback with nonzero result flow bins: same Mean/Std Dev as ROOT's ProjectionX
        # (1.625, 0.5994789404); ROOT 6.40's projection sets the in-range effective entries
        # 3.878787879, while ROOT's ResetStats of the drawn result counts the flow bins too.
        last_box = boxes[end]
        @test isapprox(last_box.values["mean1"], 1.625; rtol=1e-9)
        @test isapprox(last_box.values["std1"], 0.5994789404; rtol=1e-9)
        p = cases[end]
        all_w = [p.underflow[1]; p.sumw; p.overflow[1]]
        all_v = [p.underflow[2]; p.sumw2; p.overflow[2]]
        @test any(v -> isapprox(last_box.values["entries"], v; rtol=1e-9),
                  (sum(all_w)^2 / sum(all_v), 3.878787879))
        @test parse(Int, box_value(boxes[1], "Entries")) == 7
        @test parse(Int, box_value(boxes[2], "Entries")) == 10
    else
        @test_skip false
    end
end
