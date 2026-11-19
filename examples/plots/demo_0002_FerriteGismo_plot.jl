using ComodoIGA
using ComodoIGA.Comodo
using ComodoIGA.GeometryBasics
using ComodoIGA.Ferrite
using ComodoIGA.FerriteGismo          # re-exports TinyGismo (create..., THBSpline, RefinementBox, ...)
using ComodoIGA.GLMakie

# with adaptivity, it is better to use AT1 model

# Tuesday , 14:30 
show_control_points = true          # true / false

GLMakie.closeall()

for testCase = 1:5

    # --------------------------------------------------------------
    # Grid
    # order and nels: one number for all directions, or one value per direction
    # IGAGrid{S, R}: S = spatial dim, R = parametric dim
    # --------------------------------------------------------------
    if testCase == 1        # rectangle, different order per direction
        nels = (4, 4)
        order = (2, 2)
        geometry = linear_rectangle(0, 0, 1, 1)
        set_order_and_elements!(geometry, order, nels)
        grid = IGAGrid{2,2}(geometry)
        name = "Rectangle"

    elseif testCase == 2    # quarter annulus (arc direction needs order >= 2)
        nels = (30, 30)
        order = 2
        geometry = createNurbsQuarterAnnulus(1.0, 2.0)
        set_order_and_elements!(geometry, order, nels)
        grid = IGAGrid{2,2}(geometry)
        name = "Quarter annulus"

    elseif testCase == 3    # sphere (surface in 3D)
        nels = (20, 20)
        order = 2
        geometry = createNurbsSphere(1.0)
        set_order_and_elements!(geometry, order, nels)
        grid = IGAGrid{3,2}(geometry)
        name = "Sphere"

    elseif testCase == 4    # THB: rectangle with one corner refined locally
        # THB splines are not rational, so start from a B-spline (not a NURBS) geometry
        nels = (4, 4)
        order = 2
        geometry = linear_rectangle(0, 0, 4, 2)
        set_order_and_elements!(geometry, order, nels)

        patch = THBSpline{2}(geometry)
        # level 2, indices on level 2: coarse cells (1:2, 1:2) -> level-2 cells (1:4, 1:4)
        refineElements!(patch, RefinementBox(2, 1:4, 1:4))
        grid = IGAGrid{2,2}(patch)
        name = "THB rectangle"

    elseif testCase == 5    # cube (solid)
        nels = (4, 3, 2)
        order = 2
        geometry = createBSplineCube(1.0)
        set_order_and_elements!(geometry, order, nels)
        grid = IGAGrid{3,3}(geometry)
        name = "Cube"
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
    fig = Figure(size=(900, 700))
    title = "$name: $(length(E)) elements, $(length(Vc)) control points, order $order"

    if D == 2                               # plate: all quads
        ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y", title=title)
        meshplot!(ax, GeometryBasics.Mesh(Ve, F), color=:lightgray, strokecolor=:black,
            strokewidth=1, shading=false, transparency=false)

    elseif R == 2                           # surface in 3D: all quads
        ax = Axis3(fig[1, 1], aspect=:data, title=title)
        meshplot!(ax, GeometryBasics.Mesh(Ve, F), color=(:lightgray, 0.3), strokecolor=:black,
            strokewidth=1, shading=true, transparency=true)

    else                                    # solid: outer faces only
        ax = Axis3(fig[1, 1], aspect=:data, title=title)
        meshplot!(ax, GeometryBasics.Mesh(Ve, Fb), color=(:lightgray, 0.3), strokecolor=:black,
            strokewidth=1, shading=true, transparency=true)
    end

    # boundary edges of plates and surfaces (skipped for THB: hanging nodes break the edge matching)
    if R == 2 && !isHierarchical(grid)
        linesegments!(ax, [Ve[i] for e in Fb for i in e], color=:blue, linewidth=3)
    end

    if show_control_points
        scatter!(ax, Vc, color=:red, markersize=10, overdraw=true)
    end

    display(GLMakie.Screen(), fig)
end