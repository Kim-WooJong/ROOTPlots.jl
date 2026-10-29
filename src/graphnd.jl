"""Owned rectangular samples of z(x,y). `values[ix,iy]` belongs to the cell between adjacent edges."""
struct Grid2D
    xedges::Vector{Float64}
    yedges::Vector{Float64}
    values::Matrix{Float64}
    function Grid2D(xedges::AbstractVector{<:Real}, yedges::AbstractVector{<:Real},
                    values::AbstractMatrix{<:Real})
        grid = new(Float64.(xedges), Float64.(yedges), Matrix{Float64}(values))
        validate_grid(grid)
        return grid
    end
end

function validate_grid(grid::Grid2D)
    nx, ny = size(grid.values)
    1 <= nx <= 2048 && 1 <= ny <= 2048 && nx * ny <= 1_000_000 ||
        throw(ArgumentError("grid needs 1:1,000,000 cells with at most 2048 per axis"))
    length(grid.xedges) == nx + 1 && length(grid.yedges) == ny + 1 ||
        throw(ArgumentError("grid values must have shape (length(xedges)-1, length(yedges)-1)"))
    all(isfinite, grid.xedges) && all(isfinite, grid.yedges) &&
        all(diff(grid.xedges) .> 0) && all(diff(grid.yedges) .> 0) &&
        all(isfinite, grid.values) || throw(ArgumentError("grid needs finite, increasing edges and finite values"))
    return grid
end

"""Owned regular 3-D histogram cells. `values[ix,iy,iz]` is the content of one TH3 bin."""
struct Grid3D
    xedges::Vector{Float64}
    yedges::Vector{Float64}
    zedges::Vector{Float64}
    values::Array{Float64,3}
    function Grid3D(xedges::AbstractVector{<:Real}, yedges::AbstractVector{<:Real},
                    zedges::AbstractVector{<:Real}, values::AbstractArray{<:Real,3})
        grid = new(Float64.(xedges), Float64.(yedges), Float64.(zedges),
                   Array{Float64,3}(values))
        validate_grid(grid)
        return grid
    end
end

function validate_grid(grid::Grid3D)
    nx, ny, nz = size(grid.values)
    all(n -> 1 <= n <= 256, (nx,ny,nz)) && nx*ny*nz <= 1_000_000 ||
        throw(ArgumentError("3-D grid needs at most 256 bins per axis and 1,000,000 cells"))
    for (edges, n) in ((grid.xedges,nx),(grid.yedges,ny),(grid.zedges,nz))
        length(edges) == n+1 && all(isfinite,edges) && all(diff(edges) .> 0) ||
            throw(ArgumentError("3-D grid edges must be finite and match the cell shape"))
    end
    all(isfinite, grid.values) || throw(ArgumentError("3-D grid needs finite cell values"))
    return grid
end

"""Owned XYZ points for a genuine 3-D scatter plot; point order and duplicate xy pairs are retained."""
struct Graph3D
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    style::GraphStyle
    function Graph3D(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                     z::AbstractVector{<:Real}, style::GraphStyle)
        cloud = new(Float64.(x), Float64.(y), Float64.(z), style)
        validate_cloud(cloud)
        return cloud
    end
end

function Graph3D(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                 z::AbstractVector{<:Real}; style=GraphStyle())
    return Graph3D(x, y, z, style)
end

function validate_cloud(cloud::Graph3D)
    n = length(cloud.x)
    1 <= n <= 1_000_000 && n == length(cloud.y) == length(cloud.z) &&
        all(isfinite, cloud.x) && all(isfinite, cloud.y) && all(isfinite, cloud.z) ||
        throw(ArgumentError("XYZ needs equal, nonempty, finite coordinate arrays"))
    minimum(cloud.x) < maximum(cloud.x) && minimum(cloud.y) < maximum(cloud.y) ||
        throw(ArgumentError("XYZ rendering needs varying x and y coordinates"))
    check_graph_style(cloud.style)
    return cloud
end

"""Owned XYZ samples with symmetric absolute errors, drawn by ROOT TGraph2DErrors."""
struct Graph3DErrors
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    ex::Vector{Float64}
    ey::Vector{Float64}
    ez::Vector{Float64}
    style::GraphStyle
end

function Graph3DErrors(x::AbstractVector{<:Real},y::AbstractVector{<:Real},z::AbstractVector{<:Real};
                       xerrors=nothing,yerrors=nothing,zerrors=nothing,style=GraphStyle())
    base=Graph3D(x,y,z;style)
    n=length(base.x)
    errors(input)=input===nothing ? zeros(n) : Float64.(input)
    cloud=Graph3DErrors(base.x,base.y,base.z,errors(xerrors),errors(yerrors),errors(zerrors),style)
    validate_cloud(cloud)
    return cloud
end

function validate_cloud(cloud::Graph3DErrors)
    n=length(cloud.x)
    1<=n<=1_000_000 && all(v->length(v)==n,(cloud.y,cloud.z,cloud.ex,cloud.ey,cloud.ez)) ||
        throw(ArgumentError("XYZ error arrays must have equal nonempty lengths"))
    for (values,errors) in ((cloud.x,cloud.ex),(cloud.y,cloud.ey),(cloud.z,cloud.ez))
        all(isfinite,values) && all(isfinite,errors) && all(>=(0),errors) &&
            all(isfinite,values.-errors) && all(isfinite,values.+errors) ||
            throw(ArgumentError("XYZ coordinates and absolute errors must be finite and nonnegative"))
    end
    minimum(cloud.x)<maximum(cloud.x) && minimum(cloud.y)<maximum(cloud.y) ||
        throw(ArgumentError("XYZ rendering needs varying x and y coordinates"))
    check_graph_style(cloud.style)
    return cloud
end

"""Reject ambiguous z(x,y) samples before asking ROOT for Delaunay interpolation."""
function validate_triangulation(cloud)
    n=length(cloud.x)
    n>=3 || throw(ArgumentError("interpolation requires at least three XYZ samples"))
    length(Set(zip(cloud.x,cloud.y)))==n ||
        throw(ArgumentError("interpolation requires unique (x,y) sample coordinates"))
    # Normalize each axis before testing area to avoid overflow and scale dependence.
    xmin,xmax=extrema(cloud.x); ymin,ymax=extrema(cloud.y)
    dx=xmax-xmin; dy=ymax-ymin
    isfinite(dx) && isfinite(dy) && dx>0 && dy>0 ||
        throw(ArgumentError("interpolation coordinate spans must be finite and positive"))
    x=(cloud.x.-xmin)./dx; y=(cloud.y.-ymin)./dy
    far=argmax((x.-x[1]).^2 .+ (y.-y[1]).^2)
    area=maximum(abs.((x[far]-x[1]).*(y.-y[1]) .- (y[far]-y[1]).*(x.-x[1])))
    area>64eps(Float64) || throw(ArgumentError("interpolation requires non-collinear (x,y) samples"))
    return cloud
end

"""
    Graph2DSpec(; title="", xlabel="x", ylabel="y", zlabel="z", draw=:heatmap,
                palette=:viridis, contours=10, logx=false, logy=false, logz=false,
                size=(800, 650), xaxis=AxisStyle(), yaxis=AxisStyle(), zaxis=AxisStyle(),
                canvas=CanvasStyle())

A two-coordinate scalar-field view (colored cells, contours, or a ROOT option
string such as `"COLZ"`), used for `Grid2D`, `Hist2D`, `Profile2D`,
`PolygonHist2D`, and TF2 `ModelSurface`.

- `title`, `xlabel`, `ylabel`, `zlabel`: plot title and axis titles (`zlabel` titles
  the palette axis).
- `draw`: `:heatmap` (`COLZ`), `:contour` (`CONT3`), or a ROOT draw option string.
- `palette`: `:viridis` or `:bird` (`gStyle->SetPalette`); `contours`: number of
  contour levels (`TH1::SetContour`).
- `logx`, `logy`, `logz`: `TPad::SetLogx/SetLogy/SetLogz` (`logz` = logarithmic color scale).
- `size`: canvas size in pixels.
- `xaxis`, `yaxis`, `zaxis`: [`AxisStyle`](@ref) for the histogram axes; `zaxis` also
  styles the palette axis (ROOT's palette imports the z-axis attributes).
- `canvas`: [`CanvasStyle`](@ref) (margins override the package's layout defaults,
  grids, ROOT title box).
- `style`, `cms_label`, `atlas_label`: ROOT style and experiment label, as in
  [`GraphSpec`](@ref). With a style, flat 2-D views keep the style's margins
  (a palette still gets the package's right margin, sized for the style's z label and title).
"""
Base.@kwdef struct Graph2DSpec
    title::String = ""
    xlabel::String = "x"
    ylabel::String = "y"
    zlabel::String = "z"
    draw::Union{Symbol,String} = :heatmap
    palette::Symbol = :viridis
    contours::Int = 10
    logx::Bool = false
    logy::Bool = false
    logz::Bool = false
    size::Tuple{Int,Int} = (800,650)
    xaxis::AxisStyle = AxisStyle()
    yaxis::AxisStyle = AxisStyle()
    zaxis::AxisStyle = AxisStyle()
    canvas::CanvasStyle = CanvasStyle()
    style::Union{Nothing,Symbol} = nothing
    cms_label::Union{Nothing,Bool,NamedTuple} = nothing
    atlas_label::Union{Nothing,Bool,NamedTuple} = nothing
end

# The original eleven-field positional form remains valid.
Graph2DSpec(title,xlabel,ylabel,zlabel,draw,palette,contours,logx,logy,logz,size) =
    Graph2DSpec(;title,xlabel,ylabel,zlabel,draw,palette,contours,logx,logy,logz,size)

"""
    Graph3DSpec(; title="", xlabel="x", ylabel="y", zlabel="z", draw=:surface,
                palette=:viridis, theta=30, phi=30, logx=false, logy=false, logz=false,
                size=(900, 700), interpolate=false, resolution=(40, 40), contours=10,
                xaxis=AxisStyle(), yaxis=AxisStyle(), zaxis=AxisStyle(), canvas=CanvasStyle())

A native 3-D or scattered scalar-field view (surfaces, TH3 volumes, XYZ graphs,
TF2/TF3). Graph3D/Graph3DErrors surfaces, triangles, heatmaps and contours require
`interpolate=true`; `resolution=(nx,ny)` controls ROOT's interpolation grid.
Points/lines do not interpolate and preserve duplicate xy samples.

- `theta`, `phi`: view angles (`TPad::SetTheta/SetPhi`).
- `logx`, `logy`, `logz`: `TPad::SetLogx/SetLogy/SetLogz` (not for XYZ graphs and TF3,
  whose batch painters do not support log coordinates).
- `xaxis`, `yaxis`, `zaxis`: [`AxisStyle`](@ref) for the three axes.
- `canvas`: [`CanvasStyle`](@ref).
- `style`, `cms_label`, `atlas_label`: as in [`Graph2DSpec`](@ref) (3-D views keep
  the package's margins).
Other fields as in [`Graph2DSpec`](@ref).
"""
Base.@kwdef struct Graph3DSpec
    title::String = ""
    xlabel::String = "x"
    ylabel::String = "y"
    zlabel::String = "z"
    draw::Union{Symbol,String} = :surface
    palette::Symbol = :viridis
    theta::Float64 = 30.0
    phi::Float64 = 30.0
    logx::Bool = false
    logy::Bool = false
    logz::Bool = false
    size::Tuple{Int,Int} = (900,700)
    interpolate::Bool = false
    resolution::Tuple{Int,Int} = (40,40)
    contours::Int = 10
    xaxis::AxisStyle = AxisStyle()
    yaxis::AxisStyle = AxisStyle()
    zaxis::AxisStyle = AxisStyle()
    canvas::CanvasStyle = CanvasStyle()
    style::Union{Nothing,Symbol} = nothing
    cms_label::Union{Nothing,Bool,NamedTuple} = nothing
    atlas_label::Union{Nothing,Bool,NamedTuple} = nothing
end

Graph3DSpec(title,xlabel,ylabel,zlabel,draw,palette,theta,phi,logx,logy,logz,size) =
    Graph3DSpec(;title,xlabel,ylabel,zlabel,draw,palette,theta,phi,logx,logy,logz,size)

has_scene(spec::Union{Graph2DSpec,Graph3DSpec}) = spec.xaxis != AxisStyle() || spec.yaxis != AxisStyle() ||
    spec.zaxis != AxisStyle() || spec.canvas != CanvasStyle() || has_style_or_label(spec)

function validate_nd_spec(spec)
    all(n -> 1 <= n <= 10000, spec.size) || throw(ArgumentError("canvas size must be 1:10000"))
    validate_axis(spec.xaxis; name=:x); validate_axis(spec.yaxis; name=:y); validate_axis(spec.zaxis; name=:z)
    validate_canvas(spec.canvas)
    validate_style_and_label(spec)
    spec.palette in (:bird, :viridis) || throw(ArgumentError("palette must be :bird or :viridis"))
    for label in (spec.title, spec.xlabel, spec.ylabel, spec.zlabel)
        clean_field(label)
        occursin('\0', label) && throw(ArgumentError("NUL is not permitted in labels"))
    end
    if spec isa Graph2DSpec
        if spec.draw isa AbstractString
            root_draw_option(spec.draw)
        else
            spec.draw in (:heatmap, :contour) ||
                throw(ArgumentError("2-D draw must be :heatmap, :contour, or a ROOT option string"))
        end
        1 <= spec.contours <= 255 || throw(ArgumentError("contours must be 1:255"))
    else
        if spec.draw isa AbstractString
            root_draw_option(spec.draw)
        else
            spec.draw in (:surface,:points,:line,:triangles,:heatmap,:contour) ||
                throw(ArgumentError("invalid 3-D draw mode"))
        end
        isfinite(spec.theta) && isfinite(spec.phi) && 0 <= spec.theta <= 180 &&
            -180 <= spec.phi <= 180 || throw(ArgumentError("invalid 3-D view angles"))
        all(n->4<=n<=500,spec.resolution) && 1<=spec.contours<=255 ||
            throw(ArgumentError("XYZ interpolation resolution must be 4:500 per axis and contours 1:255"))
    end
    return spec
end

function write_nd_labels(io, spec)
    println(io, join(("LABELS", (clean_field(root_latex(label)) for label in
                 (spec.title,spec.xlabel,spec.ylabel,spec.zlabel))...), '\t'))
    write_scene(io, spec) # axis/canvas styles; nothing is written for default styles
end

function write_grid_spec(path, grid::Grid2D, spec)
    nx,ny = size(grid.values)
    mode = spec isa Graph2DSpec ? "2d" : "3d"
    raw = spec.draw isa AbstractString || spec.logx || spec.logy || spec.logz
    draw = if spec.draw isa AbstractString
        String(spec.draw)
    elseif raw
        spec.draw == :heatmap ? "COLZ" : spec.draw == :contour ? "CONT3" : "SURF2Z"
    else
        String(spec.draw)
    end
    contours = spec isa Graph2DSpec ? spec.contours : 10
    theta = spec isa Graph3DSpec ? spec.theta : 30.0
    phi = spec isa Graph3DSpec ? spec.phi : 30.0
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        if raw
            println(io,join(("DRAWGRIDRAW",spec.size...,mode,root_draw_option(draw),spec.palette,
                             contours,theta,phi,Int(spec.logx),Int(spec.logy),Int(spec.logz)),'\t'))
        else
            println(io,join(("DRAWGRID",spec.size...,mode,draw,spec.palette,contours,theta,phi),'\t'))
        end
        write_nd_labels(io,spec)
        println(io,"GRID\t",nx,'\t',ny)
        for (i,x) in enumerate(grid.xedges) println(io,"XEDGE\t",i,'\t',x) end
        for (i,y) in enumerate(grid.yedges) println(io,"YEDGE\t",i,'\t',y) end
        for iy in 1:ny, ix in 1:nx
            println(io,"CELL\t",ix,'\t',iy,'\t',grid.values[ix,iy])
        end
        println(io,"END\t1")
    end
end

function write_grid_overlay_spec(path,grids,spec,options,styles)
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        println(io,join(("DRAWGRIDOVERLAY",spec.size...,spec.palette,spec.contours,
                         Int(spec.logx),Int(spec.logy),Int(spec.logz),length(grids)),'\t'))
        write_nd_labels(io,spec)
        for (k,grid) in enumerate(grids)
            style=styles[k]
            nx,ny=size(grid.values)
            println(io,join(("GLAYER",k,options[k],style.color,style.marker_style,
                             style.marker_size,style.line_width,style.line_style),'\t'))
            println(io,"GRID\t",nx,'\t',ny)
            for (i,x) in enumerate(grid.xedges) println(io,"XEDGE\t",i,'\t',x) end
            for (i,y) in enumerate(grid.yedges) println(io,"YEDGE\t",i,'\t',y) end
            for iy in 1:ny, ix in 1:nx
                println(io,"CELL\t",ix,'\t',iy,'\t',grid.values[ix,iy])
            end
        end
        println(io,"END\t",length(grids))
    end
    return path
end

function xyz_draw_option(cloud,spec::Graph3DSpec)
    option = if spec.draw isa String
        root_draw_option(spec.draw)
    else
        get(Dict(:points=>"P0",:line=>"LINE",:surface=>"SURF2Z",:triangles=>"TRI1",
                 :heatmap=>"COLZ",:contour=>"CONT4Z"),spec.draw,"")
    end
    isempty(option) && throw(ArgumentError("XYZ requires an explicit supported draw mode"))
    cloud isa Graph3DErrors && spec.draw==:points && (option *= " ERR")
    return option
end

xyz_interpolates(option) = any(token->occursin(token,uppercase(option)),("SURF","TRI","CONT","COL"))

function validate_xyz_plot(cloud,spec::Graph3DSpec)
    validate_cloud(cloud); validate_nd_spec(spec)
    option=xyz_draw_option(cloud,spec)
    # TGraph2D has its own painter options; histogram options are not interchangeable.
    all(token->occursin(r"^(?:P0?|PCOL|LINE|ERR|TRI[W12]?|SURF[0-7]?Z?|CONT[0-5]?Z?|COLZ?)$",token),split(uppercase(option))) ||
        throw(ArgumentError("unsupported native TGraph2D option"))
    "ERR" in split(uppercase(option)) && !(cloud isa Graph3DErrors) &&
        throw(ArgumentError("ERR requires Graph3DErrors"))
    if xyz_interpolates(option)
        spec.interpolate || throw(ArgumentError("scattered interpolation requires interpolate=true"))
        validate_triangulation(cloud)
    end
    # ROOT's TGraph2D point/error painter does not implement log coordinates reliably.
    !(spec.logx || spec.logy || spec.logz) ||
        throw(ArgumentError("XYZ drawing does not support logarithmic axes"))
    check_spec_bins(spec, "an XYZ graph")
    return option
end

function write_xyz_spec(path, cloud::Union{Graph3D,Graph3DErrors}, spec::Graph3DSpec)
    option=validate_xyz_plot(cloud,spec)
    extended=cloud isa Graph3DErrors || spec.draw!=:points
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        header=("DRAWXYZ",spec.size...,cloud.style.color,cloud.style.marker_style,
                cloud.style.marker_size,spec.theta,spec.phi)
        if extended
            header=(header...,option,spec.palette,spec.resolution...,spec.contours,
                    Int(spec.logx),Int(spec.logy),Int(spec.logz),
                    cloud isa Graph3DErrors ? "symmetric" : "plain",Int(spec.interpolate),
                    cloud.style.line_width,cloud.style.line_style)
        end
        println(io,join(header,'\t'))
        write_nd_labels(io,spec)
        println(io,"XYZ\t",length(cloud.x))
        for i in eachindex(cloud.x)
            point=cloud isa Graph3DErrors ?
                ("POINT3E",i,cloud.x[i],cloud.y[i],cloud.z[i],cloud.ex[i],cloud.ey[i],cloud.ez[i]) :
                ("POINT3",i,cloud.x[i],cloud.y[i],cloud.z[i])
            println(io,join(point,'\t'))
        end
        println(io,"END\t1")
    end
end

function write_volume_spec(path, grid::Grid3D, spec::Graph3DSpec)
    nx,ny,nz = size(grid.values)
    open(path,"w") do io
        println(io,"ROOTPlotsGraph\t1")
        println(io,join(("DRAWVOLUME",spec.size...,root_draw_option(spec.draw),spec.palette,
                         spec.theta,spec.phi,Int(spec.logx),Int(spec.logy),Int(spec.logz)),'\t'))
        write_nd_labels(io,spec)
        println(io,"VOLUME\t",nx,'\t',ny,'\t',nz)
        for (i,x) in enumerate(grid.xedges) println(io,"XEDGE\t",i,'\t',x) end
        for (i,y) in enumerate(grid.yedges) println(io,"YEDGE\t",i,'\t',y) end
        for (i,z) in enumerate(grid.zedges) println(io,"ZEDGE\t",i,'\t',z) end
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            println(io,"CELL3\t",ix,'\t',iy,'\t',iz,'\t',grid.values[ix,iy,iz])
        end
        println(io,"END\t1")
    end
    return path
end

function render_nd(data, spec, outstem, paths; root_config)
    config=find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe=graphfit_executable(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        input=joinpath(dir,"graph.tsv")
        if data isa Grid2D
            write_grid_spec(input,data,spec)
        elseif data isa Grid3D
            write_volume_spec(input,data,spec)
        else
            write_xyz_spec(input,data,spec)
        end
        stem=joinpath(dir,"rendered")
        run(Cmd([exe,input,stem]))
        move_plot_outputs(stem, paths)
    end
    return paths
end

"""Draw a genuine TH3D volume with an unchanged ROOT draw-option string."""
function render_root(grid::Grid3D, spec::Graph3DSpec, outstem::AbstractString;
                     format=:png, root_config=nothing)
    paths = plot_output_paths(outstem, format)
    validate_grid(grid); validate_nd_spec(spec)
    spec.draw isa AbstractString ||
        throw(ArgumentError("Grid3D requires a ROOT option string, e.g. draw=\"BOX2Z\""))
    for (logscale, edges, axis) in ((spec.logx,grid.xedges,:x),
                                    (spec.logy,grid.yedges,:y),
                                    (spec.logz,grid.zedges,:z))
        logscale && first(edges) <= 0 &&
            throw(ArgumentError("log $axis needs positive grid edges"))
    end
    check_spec_bins(spec, "Grid3D", size(grid.values)...)
    return render_nd(grid,spec,outstem,paths;root_config)
end

"""Render Grid2D as a 2-D heatmap or contours using native ROOT."""
function render_root(grid::Grid2D, spec::Graph2DSpec, outstem::AbstractString;
                     format=:png, root_config=nothing)
    paths = plot_output_paths(outstem, format)
    validate_grid(grid); validate_nd_spec(spec)
    spec.logx && first(grid.xedges) <= 0 && throw(ArgumentError("log x needs positive grid edges"))
    spec.logy && first(grid.yedges) <= 0 && throw(ArgumentError("log y needs positive grid edges"))
    check_spec_bins(spec, "Grid2D", size(grid.values)...)
    return render_nd(grid,spec,outstem,paths;root_config)
end

"""Overlay several TH2D objects on one canvas, forwarding each ROOT option unchanged."""
function render_root(grids::AbstractVector{Grid2D}, spec::Graph2DSpec,
                     outstem::AbstractString; draw,
                     styles=fill(GraphStyle(),length(grids)), format=:png,
                     root_config=nothing)
    paths=plot_output_paths(outstem,format)
    2 <= length(grids) <= 32 || throw(ArgumentError("2-D overlay needs 2:32 layers"))
    length(draw) == length(styles) == length(grids) ||
        throw(ArgumentError("2-D overlay options and styles must match grid count"))
    validate_nd_spec(spec)
    options=String[]
    for (k,grid) in enumerate(grids)
        validate_grid(grid)
        check_graph_style(styles[k])
        spec.logx && first(grid.xedges) <= 0 &&
            throw(ArgumentError("log x needs positive grid edges"))
        spec.logy && first(grid.yedges) <= 0 &&
            throw(ArgumentError("log y needs positive grid edges"))
        option=draw[k]
        option isa AbstractString || throw(ArgumentError("2-D overlay options must be ROOT strings"))
        value=root_draw_option(option;allow_same=k>1)
        k>1 && !occursin("SAME",uppercase(value)) &&
            throw(ArgumentError("each 2-D overlay layer after the first needs SAME"))
        push!(options,value)
    end
    check_spec_bins(spec, "a 2-D overlay", size(first(grids).values)...)
    config=find_root_config(root_config)
    config === nothing && throw(ArgumentError("local native ROOT not found"))
    exe=graphfit_executable(config)
    mkpath(dirname(outstem))
    mktempdir() do dir
        input=write_grid_overlay_spec(joinpath(dir,"overlay.tsv"),grids,spec,options,styles)
        stem=joinpath(dir,"rendered")
        run(Cmd([exe,input,stem]))
        move_plot_outputs(stem,paths)
    end
    return paths
end

"""Render a gridded 3-D surface using native ROOT."""
function render_root(grid::Grid2D, spec::Graph3DSpec, outstem::AbstractString;
                     format=:png, root_config=nothing)
    paths = plot_output_paths(outstem, format)
    validate_grid(grid); validate_nd_spec(spec)
    spec.logx && first(grid.xedges) <= 0 && throw(ArgumentError("log x needs positive grid edges"))
    spec.logy && first(grid.yedges) <= 0 && throw(ArgumentError("log y needs positive grid edges"))
    spec.draw isa AbstractString || spec.draw == :surface ||
        throw(ArgumentError("Grid2D requires draw=:surface or a ROOT option string in 3-D"))
    if spec.draw == :surface
        all(n -> n >= 2, size(grid.values)) ||
            throw(ArgumentError("3-D surface needs at least two cells along both x and y"))
    end
    check_spec_bins(spec, "Grid2D", size(grid.values)...)
    return render_nd(grid,spec,outstem,paths;root_config)
end

"""Render native XYZ points/errors/lines, or explicitly opt into scattered ROOT interpolation."""
function render_root(cloud::Union{Graph3D,Graph3DErrors}, spec::Graph3DSpec, outstem::AbstractString;
                     format=:png, root_config=nothing)
    paths = plot_output_paths(outstem, format)
    validate_xyz_plot(cloud,spec)
    return render_nd(cloud,spec,outstem,paths;root_config)
end
