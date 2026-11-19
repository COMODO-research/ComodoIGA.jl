#
# functions.jl — included by module ComodoIGA
#
# Turns a FerriteIGA `BezierGrid` or a FerriteGismo `IGAGrid` into Comodo data:
#
#   E, Ve, F, Fb, Vc = IGAToComodo(grid)
#
# The grid itself is built in the demo with the packages' own calls:
#   FerriteIGA   : patch = generate_nurbs_patch(shape, nels, order; ...);  grid = BezierGrid(patch)
#                  nels is a tuple; order is ONE number (BezierGrid asserts equal orders)
#   FerriteGismo : geometry = create...(...);  set_order_and_elements!(geometry, order, nels)
#                  grid = IGAGrid{S, R}(geometry)
#                  order and nels: one number, or one value per direction
#
#   E  : elements          Hex8 (solids) or QuadFace (plates, surfaces such as a sphere)
#   Ve : element vertices  exact physical corners of the knot spans
#   F  : all faces         solids: 6 faces per element; quads: the elements themselves
#   Fb : boundary          solids: outer faces (QuadFace); quads: outer edges (LineFace)
#   Vc : control points
#
# Plot with  GeometryBasics.Mesh(Ve, F)  (all faces)  or  GeometryBasics.Mesh(Ve, Fb)  (skin of a solid).
#

using Comodo.GeometryBasics: Point, QuadFace, LineFace

const TinyGismo = FerriteGismo.TinyGismo

export IGAToComodo,
       set_order_and_elements!, linear_rectangle,
       control_points, element_mesh, mesh_faces, boundary_faces, boundary_edges,
       spatial_dim, parametric_dim, as_tuple,get_boundary_iga, get_iga_facetset,
       vertex_quadrature_rule, solution_at_vertices, displacement_at_vertices


# ==================================================================
# Helpers
# ==================================================================

"""
    as_tuple(x, n) -> NTuple{n,Int}

`x` is either one number used for all `n` directions, or one value per direction.
Used by `set_order_and_elements!` (FerriteGismo).
"""
as_tuple(x::Integer, n::Integer) = ntuple(_ -> Int(x), n)
function as_tuple(x, n::Integer)
    length(x) == n ||
        throw(ArgumentError("expected $n entries (one per parametric direction), got $(length(x))"))
    return ntuple(i -> Int(x[i]), n)
end

"""
    spatial_dim(grid) -> Int

Dimension of the space the geometry lives in (a sphere surface: 3).
"""
spatial_dim(::Ferrite.AbstractGrid{D}) where {D} = D

"""
    parametric_dim(grid) -> Int

Number of parametric directions (a sphere surface: 2, a cube: 3).
"""
parametric_dim(grid::FerriteIGA.BezierGrid) = _refdim(Ferrite.getrefshape(grid.cells[1]))
parametric_dim(::FerriteGismo.IGAGrid{S, R}) where {S, R} = R

_refdim(::Type{Ferrite.RefLine})          = 1
_refdim(::Type{Ferrite.RefQuadrilateral}) = 2
_refdim(::Type{Ferrite.RefHexahedron})    = 3

"""
    PointMerger{D}(digits)

Callable that returns one index per distinct point, so corners shared between
neighbouring elements become one vertex. Coordinates are compared after rounding
to `digits` decimals.
"""
struct PointMerger{D}
    points::Vector{Point{D,Float64}}
    lookup::Dict{NTuple{D,Float64},Int}
    digits::Int
end
PointMerger{D}(digits::Int) where {D} =
    PointMerger{D}(Point{D,Float64}[], Dict{NTuple{D,Float64},Int}(), digits)

function (m::PointMerger{D})(p) where {D}
    q   = Point{D,Float64}(Tuple(p))
    key = ntuple(d -> round(q[d]; digits = m.digits) + 0.0, D)   # + 0.0 turns -0.0 into 0.0
    return get!(m.lookup, key) do
        push!(m.points, q)
        length(m.points)
    end
end

_element_type(R) = R == 2 ? QuadFace{Int} : Hex8{Int}

function _check_pdim(R)
    R in (2, 3) || throw(ArgumentError("only parametric dimension 2 or 3 can be meshed, got $R"))
    return R
end


# ==================================================================
# FerriteIGA
# ==================================================================

"""
    element_mesh(grid::FerriteIGA.BezierGrid; digits = 10) -> (E, Ve)

A Bézier element passes through its corner Bézier points, so those are the exact
physical element corners for any geometry. They come first in each cell (Ferrite order).
The order of the grid does not matter here; note that FerriteIGA's `BezierGrid` only
accepts patches with the same order in every direction.
"""
function element_mesh(grid::FerriteIGA.BezierGrid; digits::Int = 10)
    R  = _check_pdim(parametric_dim(grid))
    D  = spatial_dim(grid)
    nc = 2^R
    merge = PointMerger{D}(digits)

    E = Vector{_element_type(R)}(undef, Ferrite.getncells(grid))
    for c in 1:Ferrite.getncells(grid)
        xb   = first(FerriteIGA.get_bezier_coordinates(grid, c))
        E[c] = _element_type(R)(ntuple(a -> merge(xb[a]), nc))
    end
    return E, merge.points
end


# ==================================================================
# FerriteGismo
# ==================================================================

"""
    linear_rectangle(x0, y0, x1, y1) -> TensorBSpline{2}

Bilinear single-element rectangle. `createBSplineRectangle` starts at degree 2,
so use this one when you need order 1.
"""
function linear_rectangle(x0, y0, x1, y1)
    kv  = TinyGismo.KnotVector([0.0, 0.0, 1.0, 1.0])
    cps = Float64[x0 y0
                  x1 y0
                  x0 y1
                  x1 y1]
    return TinyGismo.TensorBSpline{2}(TinyGismo.TensorBSplineBasis{2}(kv, kv), cps)
end

"""
    set_order_and_elements!(geometry, order, nels) -> geometry

FerriteGismo only. Raise the degree to `order`, then insert uniform knots so direction
`d` has `nels[d]` elements. `order` and `nels` are one number for all directions or one
value per direction (unlike FerriteIGA, a different order per direction is allowed). Degree first, knots second keeps maximum continuity. Knots the geometry
already has are not inserted again.
"""
function set_order_and_elements!(geometry, order, nels)
    R = Int(TinyGismo.parDim(geometry))          # G+Smo returns Int32
    S = Int(TinyGismo.targetDim(geometry))
    orders = as_tuple(order, R)
    nelems = as_tuple(nels, R)

    for dir in 1:R
        p = Int(TinyGismo.degree(geometry, dir))
        orders[dir] < p &&
            error("direction $dir already has degree $p; cannot lower it to $(orders[dir]) without changing the geometry")
        orders[dir] > p && TinyGismo.degreeElevate!(geometry, orders[dir] - p, dir)
    end

    existing = [FerriteGismo.breakpoints(FerriteGismo.IGAGrid{S, R}(geometry), dir) for dir in 1:R]
    for dir in 1:R, k in 1:(nelems[dir] - 1)
        ξ = k / nelems[dir]
        any(b -> isapprox(b, ξ; atol = 1e-12), existing[dir]) || TinyGismo.insertKnot!(geometry, ξ, dir)
    end

    got = ntuple(dir -> length(FerriteGismo.breakpoints(FerriteGismo.IGAGrid{S, R}(geometry), dir)) - 1, R)
    got == nelems ||
        @warn "requested $nelems elements but got $got: the starting knots are not a subset of the uniform ones" starting = ntuple(d -> length(existing[d]) - 1, R)
    return geometry
end

"""
    element_mesh(grid::FerriteGismo.IGAGrid; digits = 10) -> (E, Ve)

Exact physical element corners: from `exportGrid` for plates and surfaces, from the
knot-span corners mapped with `toPhysical` for solids. Works for THB grids too.
"""
function element_mesh(grid::FerriteGismo.IGAGrid; digits::Int = 10)
    R = _check_pdim(parametric_dim(grid))
    D = spatial_dim(grid)
    merge = PointMerger{D}(digits)

    if R == 2
        eg = FerriteGismo.exportGrid(grid; subdivision = 1)
        E  = [QuadFace{Int}(map(n -> merge(eg.nodes[n].x), cell.nodes)) for cell in eg.cells]
        return E, merge.points
    end

    E = Hex8{Int}[]
    for ks in grid.knotSpans
        l, u = ks.lower, ks.upper
        corners = ((l[1], l[2], l[3]), (u[1], l[2], l[3]), (u[1], u[2], l[3]), (l[1], u[2], l[3]),
                   (l[1], l[2], u[3]), (u[1], l[2], u[3]), (u[1], u[2], u[3]), (l[1], u[2], u[3]))
        push!(E, Hex8{Int}(map(ξ -> merge(FerriteGismo.toPhysical(grid, collect(ξ))), corners)))
    end
    return E, merge.points
end


# ==================================================================
# Common to both
# ==================================================================


"""
    control_points(grid) -> Vector{Point}

Spline control points. Both FerriteIGA and FerriteGismo store them as the grid nodes.
"""
function control_points(grid::Ferrite.AbstractGrid)
    D = spatial_dim(grid)
    return [Point{D,Float64}(Tuple(Ferrite.get_node_coordinate(node))) for node in grid.nodes]
end

"""
    boundary_faces(F) -> Vector{QuadFace}

Faces of a solid used by exactly one element (the outer skin).
"""
function boundary_faces(F)
    F_uni, _, counts = gunique(F; return_index = Val(true),
                                  return_counts = Val(true),
                                  sort_entries = true)
    return F_uni[isone.(counts)]
end

"""
    boundary_edges(F) -> Vector{LineFace}

Edges of a quad mesh used by exactly one quad. Collapsed edges (sphere poles) are skipped.
"""
function boundary_edges(F)
    count = Dict{Tuple{Int,Int},Int}()
    edge  = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for f in F, k in 1:4
        a, b = f[k], f[mod1(k + 1, 4)]
        a == b && continue
        key = minmax(a, b)
        count[key] = get(count, key, 0) + 1
        haskey(edge, key) || (edge[key] = (a, b))
    end
    return [LineFace{Int}(edge[key]) for (key, c) in count if c == 1]
end

"""
    mesh_faces(E) -> (F, Fb)

Solids (`Hex8`): all element faces and the outer faces.
Plates and surfaces (`QuadFace`): the quads themselves and the outer edges.
"""
function mesh_faces(E::AbstractVector{<:Hex8})
    F = element2faces(E)
    return F, boundary_faces(F)
end

mesh_faces(E::AbstractVector{<:QuadFace}) = (E, boundary_edges(E))

"""
    IGAToComodo(grid; digits = 10) -> (E, Ve, F, Fb, Vc)

Comodo data for a FerriteIGA `BezierGrid` or a FerriteGismo `IGAGrid`.

    E  : elements (Hex8 or QuadFace)
    Ve : element vertices (exact physical corners)
    F  : all faces
    Fb : boundary faces (solids) or boundary edges (plates, surfaces)
    Vc : control points
"""
function IGAToComodo(grid; digits::Int = 10)
    E, Ve = element_mesh(grid; digits)
    F, Fb = mesh_faces(E)
    Vc    = control_points(grid)
    return E, Ve, F, Fb, Vc
end

function get_boundary_iga(grid, facetset; digits = 10)
    E, Ve = element_mesh(grid; digits)
    R = parametric_dim(grid)

    refshape = R == 2 ? Ferrite.RefQuadrilateral : Ferrite.RefHexahedron
    facets = Ferrite.reference_facets(refshape)

    indices = Int[]
    for facet in facetset
        cellid, facetid = facet
        append!(indices, [E[cellid][i] for i in facets[facetid]])
    end

    return Ve[unique(indices)]
end

function get_iga_facetset(grid::FerriteGismo.IGAGrid{S,3},
                          name::AbstractString) where {S}
    # name => (parametric direction, upper side?, Ferrite face number)
    sides = Dict(
        "left"   => (1, false, 5),
        "right"  => (1, true,  3),
        "front"  => (2, false, 2),
        "back"   => (2, true,  4),
        "bottom" => (3, false, 1),
        "top"    => (3, true,  6),
    )
    haskey(sides, name) || throw(ArgumentError("Unknown side: $name"))
    dir, upper, face = sides[name]

    spans = grid.knotSpans
    bound = upper ? maximum(k.upper[dir] for k in spans) :
                    minimum(k.lower[dir] for k in spans)

    return [
        Ferrite.FacetIndex(i, face)
        for (i, k) in enumerate(spans)
        if isapprox(upper ? k.upper[dir] : k.lower[dir], bound;
                    atol = 1e-12, rtol = 0)
    ]
end 

# ==================================================================
# Evaluating a solution at the element vertices
# ==================================================================

"""
    vertex_quadrature_rule(refshape) -> QuadratureRule
    vertex_quadrature_rule(ip::Ferrite.Interpolation) -> QuadratureRule

Rule whose points are the corners of the reference element, in the same vertex order
`element_mesh` uses, so point `a` belongs to vertex `E[cellid][a]`.

Nothing is integrated with it: it only says where the field is evaluated, so the weights
carry no meaning and are set to zero. The reference shape is taken from the interpolation
(`RefQuadrilateral` for plates and surfaces, `RefHexahedron` for solids) instead of being
written out at the call site.
"""
function vertex_quadrature_rule(refshape::Type{<:Ferrite.AbstractRefShape})
    corners = Ferrite.reference_coordinates(Lagrange{refshape, 1}())
    return QuadratureRule{refshape}(zeros(length(corners)), corners)
end

vertex_quadrature_rule(ip::Ferrite.Interpolation) = vertex_quadrature_rule(Ferrite.getrefshape(ip))

#=
Each package needs its own evaluator.

FerriteIGA: `BezierCellValues` does everything; only the values are asked for, so
gradients, hessians and detJdV are switched off.

FerriteGismo: Ferrite's `reinit!` cannot be used. It goes through a geometric mapping
whose Jacobian is singular at the element corners, so it throws `det(J) is not positive`,
and for a three-dimensional knot span it has no `areaOfKnotSpan` at all. The spline basis
is therefore evaluated straight at the corners of the knot span, the same way
`FerriteGismo.toPhysical` evaluates the geometry there.

The basis is asked for its GLOBAL active functions at the point, not for the ones the cell
stores: on the boundary of a knot span G+Smo returns the actives of a neighbouring span,
in an order that does not match the cell's dofs. Global indices sidestep that, because an
`IGADofHandler` numbers a single field's dofs as `vdim * (i - 1) + d` over the global basis
function `i` (see `_collect_cell_dofs_iga!`).
=#

struct BezierVertexEvaluator{CV, G}
    cv::CV
    grid::G
end

struct SplineVertexEvaluator{R, G}
    grid::G
    vdim::Int
    refpoints::Vector{Ferrite.Vec{R, Float64}}   # corners of the reference element
    actives::Vector{Vector{Int}}                 # global basis functions, per corner
    shapevalues::Vector{Vector{Float64}}         # their values there, per corner
end

function vertex_evaluator(grid::FerriteIGA.BezierGrid, qr, ip)
    cv = FerriteIGA.BezierCellValues(
        qr, ip;
        update_gradients = false,
        update_hessians  = false,
        update_detJdV    = false
    )
    return BezierVertexEvaluator(cv, grid)
end

function vertex_evaluator(grid::FerriteGismo.IGAGrid, qr, ip)
    refpoints = collect(Ferrite.getpoints(qr))
    R         = length(eltype(refpoints))
    n         = length(refpoints)
    return SplineVertexEvaluator{R, typeof(grid)}(
        grid, Ferrite.n_components(ip), refpoints,
        [Int[] for _ in 1:n], [Float64[] for _ in 1:n]
    )
end

prepare_cell!(ev::BezierVertexEvaluator, cellid::Int) =
    Ferrite.reinit!(ev.cv, Ferrite.getcoordinates(ev.grid, cellid))

function prepare_cell!(ev::SplineVertexEvaluator{R}, cellid::Int) where {R}
    basis     = FerriteGismo.TinyGismo.basis(ev.grid.geometry)
    knotspan  = ev.grid.knotSpans[cellid]
    actives   = FerriteGismo.gsMatrix{Int32}()
    values    = FerriteGismo.gsMatrix()

    for (a, ξ) in pairs(ev.refpoints)
        x = collect(FerriteGismo.ref_to_param(ξ, knotspan))
        FerriteGismo.active!(basis, x, actives)
        FerriteGismo.TinyGismo.eval!(basis, x, values)
        ev.actives[a]     = Int.(FerriteGismo.toVector(actives))
        ev.shapevalues[a] = FerriteGismo.toVector(values)
    end
    return ev
end

vertex_value(ev::BezierVertexEvaluator, a::Int, ue) = function_value(ev.cv, a, ue)

# `ue` here is the whole solution vector: the global basis indices index into it directly.
function vertex_value(ev::SplineVertexEvaluator, a::Int, u)
    vdim = ev.vdim
    return ntuple(vdim) do d
        v = zero(eltype(u))
        @inbounds for (N, i) in zip(ev.shapevalues[a], ev.actives[a])
            v += N * u[vdim * (i - 1) + d]
        end
        v
    end
end

# The Bézier evaluator works on the dofs of one cell; the spline one indexes the whole
# vector with global basis indices.
function cell_values_of(::BezierVertexEvaluator, dh, u, cellid, dofs)
    resize!(dofs, Ferrite.ndofs_per_cell(dh, cellid))
    Ferrite.celldofs!(dofs, dh, cellid)
    return u[dofs]
end

cell_values_of(::SplineVertexEvaluator, dh, u, cellid, dofs) = u

_as_point(v::Number) = Point{1, Float64}(v)
_as_point(v)         = Point{length(v), Float64}(Tuple(v))

"""
    solution_at_vertices(grid, dh, u, ip, E, Ve) -> Vector{Point}

Value of the field interpolated by `ip` at the element vertices `Ve` of
[`IGAToComodo`](@ref): one `Point` per vertex, ready to be added to `Ve` or used as a
colour. `u` is the full solution vector and `dh` the dof handler it came from.

The coefficients in `u` belong to control points, not to vertices, so the field has to be
evaluated at the element corners; that is what this does. It works for a FerriteIGA
`BezierGrid` and for a FerriteGismo `IGAGrid`, in 2D and 3D, for a scalar or a vector
field. `dh` must carry the single field that `ip` interpolates.

`displacement_at_vertices` is the same function under the name that reads better in a
mechanics demo.
"""
function solution_at_vertices(grid, dh, u::AbstractVector, ip, E, Ve)
    qr   = vertex_quadrature_rule(ip)
    nc   = length(Ferrite.getpoints(qr))
    vdim = Ferrite.n_components(ip)
    ev   = vertex_evaluator(grid, qr, ip)

    values = zeros(Point{vdim, Float64}, length(Ve))
    dofs   = Int[]

    for cellid in 1:Ferrite.getncells(grid)
        prepare_cell!(ev, cellid)
        ue = cell_values_of(ev, dh, u, cellid, dofs)

        for a in 1:nc
            values[E[cellid][a]] = _as_point(vertex_value(ev, a, ue))
        end
    end

    return values
end

const displacement_at_vertices = solution_at_vertices
