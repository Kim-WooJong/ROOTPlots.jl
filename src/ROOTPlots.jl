module ROOTPlots

using DataFrames: AbstractDataFrame
using SpecialFunctions: gamma_inc, erf, erfc
using SHA: sha256
using Scratch: @get_scratch!

export Hist1D, Component, BinErrors, PanelSpec, PlotStyle, comparison,
       histogram, render_root, root_available, build_root, rootmap, RootMapError,
       GraphSeries, GraphSpec, GraphStyle, FitResult, fit_root, fit_curve,
       AxisStyle, CanvasStyle, PlotLine, PlotText, LegendSpec, LegendEntry,
       Grid2D, Grid3D, Graph3D, Graph2DSpec, Graph3DSpec,
       rebin, integral, effective_entries, binned_mean, binned_std,
       binned_mean_error, histstats, fit_correlation, fit_residuals, chi2_pvalue,
       histogram_points, histogram_steps, fit_pull_points,
       Hist2D, Hist3D, projection,
       HistStack, Profile1D, Profile2D, Profile3D, PolygonHist2D, profile_stats,
       ModelCurve, function_curve, Graph3DErrors, PlotPanel, PlotFigure,
       PlotPrimitive, PlotArrow, PlotBox, PlotEllipse, SecondaryAxis,
       PolarGraph, QQGraph, SplineCurve, PieChart, EfficiencyGraph, spline_points, efficiency_points, qq_points

include("histograms.jl")
include("histogramnd.jl")
include("tables.jl")
include("comparison.jl")
include("presentation.jl")
include("root_backend.jl")
include("scene.jl")
include("graphfit.jl")
include("statistics.jl")
include("histdraw.jl")
include("graphnd.jl")
include("functionplots.jl")
include("advanced_histograms.jl")
include("primitives.jl")
include("specialgraphs.jl")
include("layouts.jl")
include("parallel.jl")

end
