using ComodoIGA
using ComodoIGA.Comodo
using ComodoIGA.GeometryBasics
using ComodoIGA.Ferrite
using ComodoIGA.FerriteGismo
using ComodoIGA.GLMakie


GLMakie.closeall()
nels = (4, 3, 2)
order = 2
geometry = createBSplineCube(1.0)
set_order_and_elements!(geometry, order, nels)
grid = IGAGrid{3,3}(geometry)
name = "Cube"
# --------------------------------------------------------------
# Data
# --------------------------------------------------------------
E, Ve, F, Fb, Vc = IGAToComodo(grid)

D = spatial_dim(grid)
R = parametric_dim(grid)
# --------------------------------------------------------------
# Plot
# --------------------------------------------------------------
fig = Figure(size=(900, 700))
title = "$name: $(length(E)) elements, $(length(Vc)) control points, order $order"

ax = AxisGeom(fig[1, 1], aspect=:data, title=title)
meshplot!(ax, GeometryBasics.Mesh(Ve, Fb), color=(:lightgray, 0.3), strokecolor=:black,
    strokewidth=1, shading=true, transparency=true)

bottom = get_iga_facetset(grid, "front")
points = get_boundary_iga(grid, bottom)
scatter!(ax, points; color=:red, markersize=18)
display(GLMakie.Screen(), fig)