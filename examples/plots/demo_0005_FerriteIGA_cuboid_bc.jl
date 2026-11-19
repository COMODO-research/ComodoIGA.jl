using ComodoIGA
using ComodoIGA.Comodo
using ComodoIGA.GeometryBasics
using ComodoIGA.Ferrite
using ComodoIGA.FerriteIGA
using ComodoIGA.GLMakie

# Cuboid (rectangular prism) with boundary conditions on all faces
GLMakie.closeall()

nels = (4, 3, 2)
order = 2
Lx, Ly, Lz = 2.0, 1.5, 1.0
patch = generate_nurbs_patch(:cube, nels, order; cornerpos=(0.0, 0.0, 0.0), size=(Lx, Ly, Lz))
grid = BezierGrid(patch)
addfacetset!(grid, "left", x->x[1]≈ 0.0)  
addfacetset!(grid, "right", x->x[1]≈ Lx)
addfacetset!(grid, "bottom", x->x[3]≈ 0.0)
addfacetset!(grid, "top", x->x[3]≈ Lz)
addfacetset!(grid, "front", x->x[2]≈ 0.0)
addfacetset!(grid, "back", x->x[2]≈ Ly)

name = "Cuboid"
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

# solid: outer faces only
ax = AxisGeom(fig[1, 1], aspect=:data, title=title)
meshplot!(ax, GeometryBasics.Mesh(Ve, Fb), color=(:lightgray, 0.3), strokecolor=:black,
    strokewidth=1, shading=true, transparency=true)

points = get_boundary_iga(grid, getfacetset(grid, "front"))
scatter!(ax, points; color = :red, markersize = 18) 
display(GLMakie.Screen(), fig)
