using ROOTPlots

"""Build synthetic named columns for one control observable; not measured data."""
function synthetic_table(edges, signal, background_a, background_b, background_c, factors)
    pred = signal .+ background_a .+ background_b .+ background_c
    observed = round.(Int, pred .* factors)
    return (; low=edges[1:end-1], high=edges[2:end], data=observed,
                     signal=signal, background_a=background_a, background_b=background_b,
              background_c=background_c)
end

"""Turn synthetic named columns into a comparison; 8–12% bands are illustrative."""
function example_panel(table::NamedTuple, xlabel; relative_uncertainty=0.08)
    edges = [Float64.(table.low); Float64(table.high[end])]
    data = Hist1D(edges, Float64.(table.data), Float64.(table.data))
    components = Component[
        Component("Background C", Hist1D(edges, table.background_c, table.background_c), "#1F77B4"),
        Component("Background B", Hist1D(edges, table.background_b, table.background_b), "#9467BD"),
        Component("Background A", Hist1D(edges, table.background_a, table.background_a), "#E377C2"),
        Component("Signal", Hist1D(edges, table.signal, table.signal), "#C5001A"),
    ]
    pred = Float64.(table.signal .+ table.background_a .+ table.background_b .+ table.background_c)
    errors = BinErrors(pred .* relative_uncertainty)
    cmp = comparison(data, components;
                     data_errors=BinErrors(sqrt.(Float64.(table.data))),
                     prediction_errors=errors,
                     uncertainty_label="Illustrative unc.")
    spec = PanelSpec(xlabel; ylabel=xlabel == "Mass (GeV)" ? "Events / 20 GeV" : "Events",
                     annotation="Synthetic example")
    return cmp, spec
end

"""Return three deterministic synthetic data/model comparisons and labels; no files or physics inputs are read."""
function demo_panels()
    mass_edges = collect(20.0:20.0:420.0)
    tt_mass = Float64[150, 280, 440, 570, 640, 610, 530, 420, 330, 255,
                      200, 155, 120, 95, 75, 60, 50, 42, 30, 90]
    single_mass = Float64[12, 15, 25, 32, 35, 34, 30, 25, 19, 16,
                          13, 10, 9, 7, 6, 5, 5, 4, 4, 7]
    dibo_mass = Float64[6, 7, 10, 12, 12, 12, 11, 9, 8, 7,
                        6, 5, 4, 3, 3, 2, 2, 2, 2, 3]
    z_mass = Float64[6, 10, 16, 22, 24, 22, 18, 15, 12, 10,
                     8, 6, 5, 4, 3, 3, 2, 2, 2, 3]
    factors = Float64[1.06, .97, 1.02, .98, 1.01, 1.00, .99, 1.02, .96, 1.01,
                      .98, 1.00, .94, 1.04, 1.08, .92, 1.11, .83, .77, .96]
    mass = synthetic_table(mass_edges, tt_mass, single_mass, dibo_mass, z_mass, factors)

    jets = synthetic_table(collect(0.5:1.0:7.5),
        Float64[1650, 1830, 900, 300, 105, 24, 7],
        Float64[270, 190, 65, 22, 8, 2, 0.5],
        Float64[65, 48, 17, 5, 2, 0.5, 0.1],
        Float64[110, 70, 20, 5, 1, 0.3, 0.1],
        Float64[.99, 1.01, 1.00, .98, 1.06, 1.12, 1.14])
    bjets = synthetic_table(collect(-0.5:1.0:2.5),
        Float64[1080, 2600, 1290], Float64[240, 280, 180],
        Float64[90, 65, 30], Float64[180, 110, 45], Float64[1.04, .92, 1.01])

    panels = [example_panel(mass, "Mass (GeV)"),
              example_panel(jets, "Number of jets"; relative_uncertainty=0.10),
              example_panel(bjets, "Number of b jets"; relative_uncertainty=0.08)]
    return first.(panels), last.(panels)
end
