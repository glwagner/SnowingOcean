using Oceananigans
using Oceananigans.Units: minute, minutes, hours
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Printf
using GLMakie

grid = RectilinearGrid(size=(128, 128, 32), extent=(128, 128, 32))

g = 9.81
a = 0.8    # m
k = 60     # m
λ = 2π / k # m⁻¹
σ = sqrt(g * k) # s⁻¹
Uˢ = a^2 * σ * k # m s⁻¹

@inline ∂z_uˢ(z, t, p) = 1 / (2 * p.k) * p.Uˢ * exp(2 * p.k * z)
stokes_drift = UniformStokesDrift(∂z_uˢ = ∂z_uˢ , parameters=(; k, Uˢ))

Qᵘ = -1e-4
N² = 2e-5
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))
coriolis = FPlane(f=1e-4) # s⁻¹

wc = ZFaceField(grid)
set!(wc, 1e-3)
fill_halo_regions!(wc)
c_forcing = AdvectiveForcing(w=wc)

model = NonhydrostaticModel(; grid, coriolis, stokes_drift,
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            tracers = (:b, :c),
                            buoyancy = BuoyancyTracer(),
                            forcing = (; c=c_forcing),
                            boundary_conditions = (; u=u_bcs))

Ξ(z) = randn() * exp(z / 4)

h₀ = 5 # m

bᴺ(z) = if z < - h₀
    N² * z
else
    - N² * h₀
end

Δz = zspacings(grid, Center())

bᵢ(x, y, z) = bᴺ(z) + 1e-1 * Ξ(z) * N² * model.grid.Lz
cᵢ(x, y, z) = if z > - Δz
    rand()
else
    0
end

u★ = sqrt(abs(Qᵘ))
uᵢ(x, y, z) = u★ * 1e-1 * Ξ(z)

set!(model, u=uᵢ, w=uᵢ, b=bᵢ, c=cᵢ)

simulation = Simulation(model, Δt=10.0, stop_time=1hours)
conjure_time_step_wizard!(simulation, cfl=0.7, max_Δt=1minute)

function progress(simulation)
    u, v, w = simulation.model.velocities

    ## Print a progress message
    msg = @sprintf("i: %04d, t: %s, Δt: %s, max|u| = (%.1e, %.1e, %.1e) ms⁻¹, wall time: %s\n",
                   iteration(simulation),
                   prettytime(time(simulation)),
                   prettytime(simulation.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   prettytime(simulation.run_wall_time))

    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))

Nz = size(grid, 3)
outputs = merge(model.velocities, model.tracers)
simulation.output_writers[:xy] = JLD2OutputWriter(model, outputs,
                                                  indices = (:, :, Nz),
                                                  schedule = TimeInterval(10minutes),
                                                  filename = "langmuir_particles_xy.jld2",
                                                  overwrite_existing = true)

simulation.output_writers[:yz] = JLD2OutputWriter(model, outputs,
                                                  indices = (1, :, :),
                                                  schedule = TimeInterval(10minutes),
                                                  filename = "langmuir_particles_yz.jld2",
                                                  overwrite_existing = true)

run!(simulation)

wt = FieldTimeSeries("langmuir_particles_xy.jld2", "w")
ct = FieldTimeSeries("langmuir_particles_xy.jld2", "c")
Nt = length(wt)

fig = Figure(size=(1200, 600))
axw = Axis(fig[1, 1])
axc = Axis(fig[1, 2])
slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value
wn = @lift interior(wt[$n], :, :, 1)
cn = @lift interior(ct[$n], :, :, 1)

wlim = 1e-2
clim = 1e1

heatmap!(axw, wn, colorrange=(-wlim, wlim), colormap=:balance)
heatmap!(axc, cn, colorrange=(0, clim), colormap=:solar)
display(fig)

