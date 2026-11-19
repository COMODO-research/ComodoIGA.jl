using ComodoIGA
using ComodoIGA.Comodo
using ComodoIGA.GeometryBasics
using ComodoIGA.Ferrite
using ComodoIGA.FerriteIGA
using ComodoIGA.GLMakie

Vec = FerriteIGA.Vec

GLMakie.closeall()

# Material: plane strain
function get_material(; E, ν)
    λ = E * ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2 * (1 + ν))

    δ(i, j) = i == j ? 1.0 : 0.0

    g(i, j, k, l) =
        λ * δ(i, j) * δ(k, l) +
        μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))

    return SymmetricTensor{4,2}(g)
end

# Local stiffness matrix
function integrate_element!(ke, C::SymmetricTensor{4,2}, cv)
    n = getnbasefunctions(cv)
    δɛ = [zero(SymmetricTensor{2,2,Float64}) for _ in 1:n]

    for qp in 1:getnquadpoints(cv)
        for i in 1:n
            δɛ[i] = symmetric(shape_gradient(cv, qp, i))
        end

        dΩ = getdetJdV(cv, qp)

        for i in 1:n, j in 1:n
            ke[i, j] += (δɛ[i] ⊡ C ⊡ δɛ[j]) * dΩ
        end
    end

    return ke
end

# Local traction force
function integrate_traction_force!(fe, traction::Vec{2}, fv)
    for qp in 1:getnquadpoints(fv)
        dA = getdetJdV(fv, qp)

        for i in 1:getnbasefunctions(fv)
            δu = shape_value(fv, qp, i)
            fe[i] += (traction ⋅ δu) * dA
        end
    end

    return fe
end

# Global assembly
function assemble_problem(dh, grid, cv, fv, C, traction)
    f = zeros(ndofs(dh))
    K = allocate_matrix(dh)
    assembler = start_assemble(K, f)

    n = getnbasefunctions(cv)
    dofs = zeros(Int, n)
    fe = zeros(n)
    ke = zeros(n, n)

    for cellid in 1:getncells(grid)
        fill!(fe, 0.0)
        fill!(ke, 0.0)

        celldofs!(dofs, dh, cellid)
        reinit!(cv, getcoordinates(grid, cellid))

        integrate_element!(ke, C, cv)
        assemble!(assembler, dofs, ke, fe)
    end

    for (cellid, faceid) in getfacetset(grid, "left")
        fill!(fe, 0.0)

        celldofs!(dofs, dh, cellid)
        reinit!(fv, getcoordinates(grid, cellid), faceid)

        integrate_traction_force!(fe, traction, fv)
        f[dofs] += fe
    end

    return K, f
end

# Grid
# order = 2
# nels = (20, 10)

order = 2
nels = (500, 500) 

color_max = 0.5 
nurbsmesh = generate_nurbs_patch(:plate_with_hole, nels, order)
grid = BezierGrid(nurbsmesh)

E, Ve, F, Fb, Vc = IGAToComodo(grid)

# Boundary sets
addfacetset!(grid, "left", x -> x[1] ≈ -4.0)
addfacetset!(grid, "bot", x -> x[2] ≈ 0.0)
addfacetset!(grid, "right", x -> x[1] ≈ 0.0)

# IGA interpolation and quadrature
ip_geo = IGAInterpolation{RefQuadrilateral,order}()
ip_u = ip_geo^2

qr_cell = QuadratureRule{RefQuadrilateral}(4)
qr_face = FacetQuadratureRule{RefQuadrilateral}(3)

cv = BezierCellValues(qr_cell, ip_u)
fv = BezierFacetValues(qr_face, ip_u)

# DOF handler
dh = DofHandler(grid)
add!(dh, :u, ip_u)
close!(dh)

# Dirichlet boundary conditions
ch = ConstraintHandler(dh)

dbc1 = Dirichlet(:u, getfacetset(grid, "bot"), (x, t) -> 0.0, 2)
dbc2 = Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> 0.0, 1)

add!(ch, dbc1)
add!(ch, dbc2)
close!(ch)
FerriteIGA.update!(ch, 0.0)

# Material and traction
stiffmat = get_material(E = 100.0, ν = 0.3)
traction = Vec((-10.0, 0.0))

# Solve
@time K, f = assemble_problem(dh, grid, cv, fv, stiffmat, traction)
apply!(K, f, ch)

@time U = K \ f

# # Displacement coefficients associated with control points
# u_control = reshape(U, 2, :)

# ux = u_control[1, :]
# uy = u_control[2, :]

# disp_control = [Point{2,Float64}(ux[j], uy[j]) for j in eachindex(ux)]

# # Physical displacement at Ve
# disp_points = displacement_at_vertices(grid, dh, U, ip_u, E, Ve)

# # Step 0 is undeformed, followed by 20 loading increments
# numSteps = 20
# scale = 1.0
# load_factors = range(0.0, 1.0; length = numSteps + 1)

# UT = [α .* disp_points for α in load_factors]
# UT_mag = [[sqrt(p[1]^2 + p[2]^2) for p in displacement] for displacement in UT]

# ut_mag_max = maximum.(UT_mag)

# # Deformed element vertices
# VT = [Ve .+ scale .* displacement for displacement in UT]

# # Deformed control points
# VcT = [Vc .+ (scale * α) .* disp_control for α in load_factors]

# # Plot limits
# xmin = minimum(p[1] for vertices in VT for p in vertices)
# xmax = maximum(p[1] for vertices in VT for p in vertices)
# ymin = minimum(p[2] for vertices in VT for p in vertices)
# ymax = maximum(p[2] for vertices in VT for p in vertices)

# padding = 0.05 * max(xmax - xmin, ymax - ymin)
# # color_max = max(maximum(ut_mag_max), eps(Float64))
# color_max = 0.5
# # Visualization setup
# fig_disp = Figure(size = (900, 700))
# stepStart = 0

# ax3 = Axis(fig_disp[1, 1], aspect = DataAspect(), xlabel = "X", ylabel = "Y", title = "Step: $stepStart")

# xlims!(ax3, xmin - padding, xmax + padding)
# ylims!(ax3, ymin - padding, ymax + padding)

# hp = meshplot!(ax3, F, VT[stepStart + 1]; 
# strokecolor = :black, strokewidth = 2,
#  shading = false, color = UT_mag[stepStart + 1], 
#  transparency = false, colormap = :turbo, colorrange = (0.0, color_max))

# Colorbar(fig_disp[1, 2], hp.plots[1], label = "Displacement magnitude")

# hSlider = Slider(fig_disp[2, 1], range = 0:numSteps, startvalue = stepStart, linewidth = 30)

# on(hSlider.value) do stepIndex
#     stateIndex = stepIndex + 1
#     hp[1] = GeometryBasics.Mesh(VT[stateIndex], F)
#     hp.color = UT_mag[stateIndex]
#     ax3.title = "Step: $stepIndex"
# end

# slidercontrol(hSlider, ax3)
# display(GLMakie.Screen(), fig_disp)

# # Export the full solution
# VTKIGAFile("plate_with_hole.vtu", grid) do vtk
#     write_solution(vtk, dh, U)
# end;