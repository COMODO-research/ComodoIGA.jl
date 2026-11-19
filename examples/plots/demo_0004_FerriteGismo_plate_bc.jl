using ComodoIGA
using ComodoIGA.Comodo
using ComodoIGA.GeometryBasics
using ComodoIGA.Ferrite
using ComodoIGA.FerriteGismo
using ComodoIGA.GLMakie

show_control_points = false          # true / false

GLMakie.closeall()

nels = (4, 4)
order = (2, 2)
geometry = linear_rectangle(0, 0, 1, 1)
set_order_and_elements!(geometry, order, nels)
grid = IGAGrid{2,2}(geometry)
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
meshplot!(ax, GeometryBasics.Mesh(Ve, F), color = (:lightgray, 0.15), strokecolor=:black, strokewidth=1, shading=false, transparency=false)

if show_control_points
    scatter!(ax, Vc, color=:red, markersize=10, overdraw= true)
end
points = get_boundary_iga(grid, getfacetset(grid, "top"))
scatter!(ax, points; color = :red, markersize = 18)  
display(GLMakie.Screen(), fig)