# Native TF2/TF3 objects built from predefined formulas, and Julia callables
# sampled on grids. Arbitrary C++/TFormula expression strings are not accepted.

const ROOT_TF2_NPAR = Dict(:gaus2d=>5, :bigaus=>6, :pol1_2d=>3, :pol2_2d=>6)
const ROOT_TF3_NPAR = Dict(:gaus3d=>7, :pol1_3d=>4, :pol2_3d=>10)

function validate_surface_parameters(model, p)
    npar = haskey(ROOT_TF2_NPAR, model) ? ROOT_TF2_NPAR[model] :
           haskey(ROOT_TF3_NPAR, model) ? ROOT_TF3_NPAR[model] :
           throw(ArgumentError("unsupported analytic model $model; choose from " *
                 "$(sort(collect(keys(ROOT_TF2_NPAR)))) or $(sort(collect(keys(ROOT_TF3_NPAR))))"))
    length(p) == npar && all(isfinite, p) ||
        throw(ArgumentError("model $model needs $npar finite parameters"))
    model in (:gaus2d, :bigaus) && !(p[3] > 0 && p[5] > 0) &&
        throw(ArgumentError("Gaussian widths must be positive"))
    model == :bigaus && !(abs(p[6]) < 1) && throw(ArgumentError("bigaus correlation must satisfy |rho| < 1"))
    model == :gaus3d && !(p[3] > 0 && p[5] > 0 && p[7] > 0) &&
        throw(ArgumentError("Gaussian widths must be positive"))
    return p
end

"""Evaluate a predefined TF2/TF3 model in Julia (used for ranges and validation only)."""
function analytic_value(model::Symbol, p, x, y, z=0.0)
    model == :gaus2d && return p[1]*exp(-0.5*((x-p[2])/p[3])^2 - 0.5*((y-p[4])/p[5])^2)
    if model == :bigaus
        u, v, r = (x-p[2])/p[3], (y-p[4])/p[5], p[6]
        return p[1]*exp(-0.5/(1-r*r)*(u*u + v*v - 2r*u*v))
    end
    model == :pol1_2d && return p[1] + p[2]*x + p[3]*y
    model == :pol2_2d && return p[1] + p[2]*x + p[3]*y + p[4]*x*x + p[5]*x*y + p[6]*y*y
    model == :gaus3d && return p[1]*exp(-0.5*((x-p[2])/p[3])^2 - 0.5*((y-p[4])/p[5])^2 - 0.5*((z-p[6])/p[7])^2)
    model == :pol1_3d && return p[1] + p[2]*x + p[3]*y + p[4]*z
    model == :pol2_3d && return p[1] + p[2]*x + p[3]*y + p[4]*z + p[5]*x*x + p[6]*y*y + p[7]*z*z +
                                p[8]*x*y + p[9]*x*z + p[10]*y*z
    throw(ArgumentError("unsupported analytic model $model"))
end

function finite_range(r, name)
    t = Tuple(Float64.(r))
    length(t) == 2 && all(isfinite, t) && t[1] < t[2] ||
        throw(ArgumentError("$name must be two finite increasing values"))
    return t
end

"""
    ModelSurface(model, parameters; xrange, yrange, npx=60, npy=60, label="",
                 style=GraphStyle(color="#000000", line_width=1))

A genuine ROOT **TF2** built from a predefined formula (the 2-D analogue of
[`ModelCurve`](@ref)); no arbitrary C++ or TFormula string is evaluated.

| `model` | parameters | function |
|---|---|---|
| `:gaus2d` | `[A, mx, sx, my, sy]` | `A exp(-((x-mx)/sx)²/2 - ((y-my)/sy)²/2)` |
| `:bigaus` | `[A, mx, sx, my, sy, rho]` | correlated Gaussian, `abs(rho) < 1` |
| `:pol1_2d` | `[p0, px, py]` | `p0 + px x + py y` |
| `:pol2_2d` | `[p0, px, py, pxx, pxy, pyy]` | full quadratic |

`npx`/`npy` set ROOT's drawing resolution (`TF2::SetNpx/SetNpy`). Draw with
`Graph2DSpec` (`draw=:heatmap` → `COLZ`, `:contour` → `CONT3`, or a TF2-compatible
ROOT option such as `"CONT1Z"`, `"COL"`) or `Graph3DSpec` (`:surface` → `SURF2Z`,
or e.g. `"SURF1"`, `"LEGO2"`). ROOT paints a TF2 by evaluating it on that grid,
so surfaces and contours are grid renderings of the analytic function.
`style.color`/`line_width`/`line_style` set the TF2 line attributes (used by
mesh lines and `CONT3`).
"""
struct ModelSurface
    model::Symbol
    parameters::Vector{Float64}
    xrange::Tuple{Float64,Float64}
    yrange::Tuple{Float64,Float64}
    npx::Int
    npy::Int
    label::String
    style::GraphStyle
end

function ModelSurface(model::Symbol, parameters::AbstractVector{<:Real}; xrange, yrange,
                      npx::Integer=60, npy::Integer=60, label="", style=GraphStyle(color="#000000", line_width=1))
    s = ModelSurface(model, Float64.(parameters), finite_range(xrange, "xrange"),
                     finite_range(yrange, "yrange"), Int(npx), Int(npy), String(label), style)
    return validate_analytic(s)
end

"""
    ModelVolume(model, parameters; xrange, yrange, zrange, level=0.0,
                npx=15, npy=15, npz=15, style=GraphStyle(color="#4477aa"))

A genuine ROOT **TF3** built from a predefined formula; no arbitrary C++ or
TFormula string is evaluated.

| `model` | parameters | function |
|---|---|---|
| `:gaus3d` | `[A, mx, sx, my, sy, mz, sz]` | 3-D Gaussian |
| `:pol1_3d` | `[p0, px, py, pz]` | `p0 + px x + py y + pz z` |
| `:pol2_3d` | `[p0, px, py, pz, pxx, pyy, pzz, pxy, pxz, pyz]` | general quadric |

Only drawings that ROOT can produce in batch mode without OpenGL are offered,
selected with `Graph3DSpec(draw=...)`:

- `:iso` — ROOT's native TF3 painter draws the implicit surface
  `f(x, y, z) = level` (internally `f - level = 0`), shaded and meshed in
  `style.color`. The implicit-surface mesh uses ROOT's own internal grid.
- `:box` — the TF3 is evaluated natively at the bin centres of an
  `npx × npy × npz` TH3D and drawn with `BOX2Z` (box size and color follow the
  value); bins with `f ≤ level` are not drawn, so both modes show `f > level`.
  This is a sampled display, not a volume rendering.

Log axes, OpenGL options (`GLISO`, `GLBOX`, ...), transparent volume rendering,
and several iso levels in one pad are not provided.
"""
struct ModelVolume
    model::Symbol
    parameters::Vector{Float64}
    xrange::Tuple{Float64,Float64}
    yrange::Tuple{Float64,Float64}
    zrange::Tuple{Float64,Float64}
    level::Float64
    npx::Int
    npy::Int
    npz::Int
    style::GraphStyle
end

function ModelVolume(model::Symbol, parameters::AbstractVector{<:Real}; xrange, yrange, zrange,
                     level::Real=0.0, npx::Integer=15, npy::Integer=15, npz::Integer=15,
                     style=GraphStyle(color="#4477aa"))
    v = ModelVolume(model, Float64.(parameters), finite_range(xrange, "xrange"), finite_range(yrange, "yrange"),
                    finite_range(zrange, "zrange"), Float64(level), Int(npx), Int(npy), Int(npz), style)
    return validate_analytic(v)
end

function validate_analytic(s::ModelSurface)
    haskey(ROOT_TF2_NPAR, s.model) || throw(ArgumentError("ModelSurface model must be one of $(sort(collect(keys(ROOT_TF2_NPAR))))"))
    validate_surface_parameters(s.model, s.parameters)
    finite_range(s.xrange, "xrange"); finite_range(s.yrange, "yrange")
    4 <= s.npx <= 1000 && 4 <= s.npy <= 1000 || throw(ArgumentError("TF2 npx/npy must be in 4:1000"))
    scene_text(s.label); check_graph_style(s.style)
    return s
end

function validate_analytic(v::ModelVolume)
    haskey(ROOT_TF3_NPAR, v.model) || throw(ArgumentError("ModelVolume model must be one of $(sort(collect(keys(ROOT_TF3_NPAR))))"))
    validate_surface_parameters(v.model, v.parameters)
    finite_range(v.xrange, "xrange"); finite_range(v.yrange, "yrange"); finite_range(v.zrange, "zrange")
    isfinite(v.level) || throw(ArgumentError("iso level must be finite"))
    all(n -> 2 <= n <= 100, (v.npx, v.npy, v.npz)) || throw(ArgumentError("TF3 sampling npx/npy/npz must be in 2:100"))
    check_graph_style(v.style)
    return v
end

const TF2_OPTION = r"^(?:COL[0-2]?Z?|CONT[0-4]?Z?|SURF[1-7]?Z?|LEGO[1-4]?Z?|TEXT|Z)$"

"""Native TF2 option for a spec; presets match Grid2D, strings are restricted to TF2 painters."""
function tf2_option(spec)
    if spec.draw isa AbstractString
        option = root_draw_option(spec.draw)
        tokens = split(uppercase(option))
        !isempty(tokens) && all(t -> occursin(TF2_OPTION, t), tokens) ||
            throw(ArgumentError("unsupported TF2 option; use COL/COLZ, CONT0-4[Z], SURF[1-7][Z], LEGO[1-4][Z]"))
        any(t -> startswith(t, "SURF") || startswith(t, "LEGO"), tokens) == (spec isa Graph3DSpec) ||
            throw(ArgumentError("SURF/LEGO need Graph3DSpec; COL/CONT need Graph2DSpec"))
        return option
    end
    spec isa Graph2DSpec && spec.draw == :heatmap && return "COLZ"
    spec isa Graph2DSpec && spec.draw == :contour && return "CONT3"
    spec isa Graph3DSpec && spec.draw == :surface && return "SURF2Z"
    throw(ArgumentError("ModelSurface supports Graph2DSpec(:heatmap/:contour), Graph3DSpec(:surface), or TF2 option strings"))
end

function tf3_mode(spec)
    spec isa Graph3DSpec && spec.draw in (:iso, :box) ||
        throw(ArgumentError("ModelVolume needs Graph3DSpec(draw=:iso) or Graph3DSpec(draw=:box)"))
    return spec.draw
end

function validate_analytic_spec(s::ModelSurface, spec)
    spec isa Union{Graph2DSpec,Graph3DSpec} || throw(ArgumentError("ModelSurface needs Graph2DSpec or Graph3DSpec"))
    validate_nd_spec(spec)
    option = tf2_option(spec)
    spec.logx && s.xrange[1] <= 0 && throw(ArgumentError("log x needs a positive x range"))
    spec.logy && s.yrange[1] <= 0 && throw(ArgumentError("log y needs a positive y range"))
    xs = range(s.xrange...; length=s.npx); ys = range(s.yrange...; length=s.npy)
    values = [analytic_value(s.model, s.parameters, x, y) for x in xs, y in ys]
    all(isfinite, values) || throw(ArgumentError("TF2 evaluates to non-finite values on its grid"))
    spec.logz && !all(>(0), values) && throw(ArgumentError("log z needs positive TF2 values on the drawing grid"))
    return option, extrema(values)
end

function validate_analytic_spec(v::ModelVolume, spec)
    spec isa Graph3DSpec || throw(ArgumentError("ModelVolume needs Graph3DSpec"))
    mode = tf3_mode(spec)
    validate_nd_spec(Graph3DSpec(; (k => getfield(spec, k) for k in fieldnames(Graph3DSpec) if k != :draw)...,
                                 draw=:surface))
    (spec.logx || spec.logy || spec.logz) && throw(ArgumentError("ModelVolume does not support log axes"))
    centers(r, n) = [r[1] + (i - 0.5) * (r[2] - r[1]) / n for i in 1:n]
    for x in centers(v.xrange, v.npx), y in centers(v.yrange, v.npy), z in centers(v.zrange, v.npz)
        isfinite(analytic_value(v.model, v.parameters, x, y, z)) ||
            throw(ArgumentError("TF3 evaluates to non-finite values on its sampling grid"))
    end
    return mode
end

function write_analytic_spec(path, s::ModelSurface, spec)
    option, (zlo, zhi) = validate_analytic_spec(s, spec)
    theta = spec isa Graph3DSpec ? spec.theta : 30.0
    phi = spec isa Graph3DSpec ? spec.phi : 30.0
    contours = spec.contours
    st = s.style
    open(path, "w") do io
        println(io, "ROOTPlotsGraph\t1")
        println(io, join(("DRAWTF2", spec.size..., option, spec.palette, contours, theta, phi,
                          Int(spec.logx), Int(spec.logy), Int(spec.logz)), '\t'))
        write_nd_labels(io, spec)
        println(io, join(("TF2", s.model, s.xrange..., s.yrange..., s.npx, s.npy, length(s.parameters),
                          st.color, st.line_width, st.line_style, s.parameters...), '\t'))
        println(io, "END\t1")
    end
    return path
end

function write_analytic_spec(path, v::ModelVolume, spec)
    mode = validate_analytic_spec(v, spec)
    st = v.style
    open(path, "w") do io
        println(io, "ROOTPlotsGraph\t1")
        println(io, join(("DRAWTF3", spec.size..., mode, spec.palette, spec.theta, spec.phi), '\t'))
        write_nd_labels(io, spec)
        println(io, join(("TF3", v.model, v.xrange..., v.yrange..., v.zrange..., v.npx, v.npy, v.npz,
                          v.level, length(v.parameters), st.color, st.line_width, st.line_style,
                          v.parameters...), '\t'))
        println(io, "END\t1")
    end
    return path
end

render_root(s::Union{ModelSurface,ModelVolume}, spec::Union{Graph2DSpec,Graph3DSpec}, outstem::AbstractString;
            format=:png, root_config=nothing) =
    render_root(PlotPanel(s, spec), outstem; format, root_config)

"""
    function_grid(f, xrange, yrange; n=(60, 60)) -> Grid2D
    function_grid(f, xrange, yrange, zrange; n=(20, 20, 20)) -> Grid3D

Sample a Julia callable `f(x, y)` or `f(x, y, z)` at the **cell centres** of a
regular grid and return a `Grid2D`/`Grid3D` (drawn as TH2D/TH3D). This is the
grid analogue of [`function_curve`](@ref): a display-only sampling of arbitrary
Julia code, not a ROOT TF2/TF3. Values must be finite reals.
"""
function function_grid(f, xrange, yrange; n=(60, 60))
    xr, yr = finite_range(xrange, "xrange"), finite_range(yrange, "yrange")
    length(n) == 2 && all(k -> k isa Integer && 1 <= k <= 2048, n) && prod(n) <= 1_000_000 ||
        throw(ArgumentError("2-D function grid needs n=(nx, ny) with 1:2048 cells per axis"))
    xe = collect(range(xr...; length=n[1]+1)); ye = collect(range(yr...; length=n[2]+1))
    xc = (xe[1:end-1] .+ xe[2:end]) ./ 2; yc = (ye[1:end-1] .+ ye[2:end]) ./ 2
    raw = [f(x, y) for x in xc, y in yc]
    all(v -> v isa Real && isfinite(v), raw) || throw(ArgumentError("function grid must return finite real values"))
    return Grid2D(xe, ye, Float64.(raw))
end

function function_grid(f, xrange, yrange, zrange; n=(20, 20, 20))
    rs = (finite_range(xrange, "xrange"), finite_range(yrange, "yrange"), finite_range(zrange, "zrange"))
    length(n) == 3 && all(k -> k isa Integer && 1 <= k <= 256, n) && prod(n) <= 1_000_000 ||
        throw(ArgumentError("3-D function grid needs n=(nx, ny, nz) with 1:256 cells per axis"))
    edges = [collect(range(r...; length=k+1)) for (r, k) in zip(rs, n)]
    c = [(e[1:end-1] .+ e[2:end]) ./ 2 for e in edges]
    raw = [f(x, y, z) for x in c[1], y in c[2], z in c[3]]
    all(v -> v isa Real && isfinite(v), raw) || throw(ArgumentError("function grid must return finite real values"))
    return Grid3D(edges..., Float64.(raw))
end
