using ComodoIGA
using ComodoIGA.Comodo
using ComodoIGA.GeometryBasics
using ComodoIGA.Ferrite
using ComodoIGA.FerriteIGA
using ComodoIGA.GLMakie

# some hints in IGA:
# grid.cells: - stores the elements.
# Each cell references the control-point indices whose basis functions influence that element.
# grid.nodes: stores the spline’s control-point coordinates.
# an element is a region between consecutive knot values—a knot span
# The associated control points, their basis functions, and any NURBS weights together define the geometry within that region.

show_control_points = true          # true / false

GLMakie.closeall()

for testCase = 1:5

    # --------------------------------------------------------------
    # Grid
    # nels : one value per direction (tuple)
    # order: one number — BezierGrid requires the same order in every direction
    # --------------------------------------------------------------
    if testCase == 1        # rectangle
        nels  = (10, 5)
        order = 2
        patch = generate_nurbs_patch(:rectangle, nels, order; size = (4.0, 2.0))
        grid  = BezierGrid(patch)
        name  = "Rectangle"

    elseif testCase == 2    # plate with hole (FerriteIGA only allows order 2)
        nels  = (10, 5)
        order = 2
        patch = generate_nurbs_patch(:plate_with_hole, nels, order)
        grid  = BezierGrid(patch)
        name  = "Plate with hole"

    elseif testCase == 3    # ring (FerriteIGA only builds it as 4×1 elements, order 2)
        nels  = (4, 1)
        order = 2
        patch = generate_nurbs_patch(:ring, nels, order; ri = 1.0, ro = 2.0)

        # refine afterwards: each pass inserts a knot in the middle of every span
        nrefine = 1
        for pass in 1:nrefine, dir in 1:2
            kv = unique(patch.knot_vectors[dir])
            for ξ in (kv[1:end-1] .+ kv[2:end]) ./ 2
                patch = knotinsertion(patch, ξ; dir)
            end
        end

        grid  = BezierGrid(patch)
        name  = "Ring"

    elseif testCase == 4    # spherical patch (surface in 3D)
        nels  = (8, 8)
        order = 2
        patch = generate_nurbs_patch(:doubly_curved, nels, order; r1 = 1.0, r2 = 1.0, α1 = π/2, α2 = π/2)
        grid  = BezierGrid(patch)
        name  = "Spherical patch"

    elseif testCase == 5    # cube (solid)
        nels  = (4, 3, 2)
        order = 2
        patch = generate_nurbs_patch(:cube, nels, order; cornerpos = (0.0, 0.0, 0.0), size = (2.0, 1.5, 1.0))
        grid  = BezierGrid(patch)
        name  = "Cube"
    end

    # --------------------------------------------------------------
    # Data
    # --------------------------------------------------------------
    E, Ve, F, Fb, Vc = IGAToComodo(grid)

    D = spatial_dim(grid)
    R = parametric_dim(grid)
    println("testCase $testCase ($name): ", length(E), " elements, ", length(Ve),
            " element vertices, ", length(Vc), " control points")

    # --------------------------------------------------------------
    # Plot
    # --------------------------------------------------------------
    fig   = Figure(size = (900, 700))
    title = "$name: $(length(E)) elements, $(length(Vc)) control points, order $order"

    if D == 2                               # plate: all quads, boundary edges in blue
        ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "X", ylabel = "Y", title = title)
        meshplot!(ax, GeometryBasics.Mesh(Ve, F), color = :lightgray, strokecolor = :black,
                  strokewidth = 1, shading = false, transparency = false)
        linesegments!(ax, [Ve[i] for e in Fb for i in e], color = :blue, linewidth = 3)

    elseif R == 2                           # surface in 3D: all quads
        ax = Axis3(fig[1, 1], aspect = :data, title = title)
        meshplot!(ax, GeometryBasics.Mesh(Ve, F), color = (:lightgray, 0.3), strokecolor = :black,
                  strokewidth = 1, shading = true, transparency = true)
        linesegments!(ax, [Ve[i] for e in Fb for i in e], color = :blue, linewidth = 3)

    else                                    # solid: outer faces only
        ax = Axis3(fig[1, 1], aspect = :data, title = title)
        meshplot!(ax, GeometryBasics.Mesh(Ve, Fb), color = (:lightgray, 0.3), strokecolor = :black,
                  strokewidth = 1, shading = true, transparency = true)
    end

    if show_control_points
        scatter!(ax, Vc, color = :red, markersize = 10, overdraw = true)
    end

    display(GLMakie.Screen(), fig)
end