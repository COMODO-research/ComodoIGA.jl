using ComodoIGA
using ComodoIGA.Comodo
using ComodoIGA.GeometryBasics
using ComodoIGA.Ferrite
using ComodoIGA.FerriteIGA
using ComodoIGA.GLMakie

show_control_points = false          # true / false

GLMakie.closeall()
nels = (10, 5)
order = 2
Lx = 4.0
Ly = 2.0
patch = generate_nurbs_patch(:rectangle, nels, order; size=(4.0, 2.0))
grid = BezierGrid(patch)
addfacetset!(grid, "left", x->x[1]≈ 0.0) 
addfacetset!(grid, "right", x->x[1]≈ Lx)  
addfacetset!(grid, "bottom", x->x[2]≈ 0.0)   
addfacetset!(grid, "top", x->x[2]≈ Ly)   
name = "Rectangle"
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
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y", title=title)
meshplot!(ax, GeometryBasics.Mesh(Ve, F), color=:lightgray, strokecolor=:black,
    strokewidth=1, shading=false, transparency=false)
# linesegments!(ax, [Ve[i] for e in Fb for i in e], color=:blue, linewidth=3)

if show_control_points
    scatter!(ax, Vc, color=:red, markersize=10, overdraw=true)
end

points = get_boundary_iga(grid, getfacetset(grid, "bottom"))
scatter!(ax, points; color = :red, markersize = 18) 
display(GLMakie.Screen(), fig)