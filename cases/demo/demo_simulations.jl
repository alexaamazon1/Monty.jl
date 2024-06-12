using DrWatson
@quickactivate "cases"
using Monty
using Random
using Distributions
using Meshes
using GeoStatsFunctions
using GeoStatsProcesses

##

import CairoMakie as Mke

##

const rng = Xoshiro(1)

##

# base grid
grid = CartesianGrid((8, 8), (0.0, 0.0), (10.0, 10.0))

# split treatment and control cells, alternating
treatment = reduce(vcat, [grid[i, :] for i ∈ 1:2:8])
control = reduce(vcat, [grid[i, :] for i ∈ 2:2:8])

# take cell centroids for each group
treatment_points = treatment .|> centroid
control_points = control .|> centroid

##

# sampling before and after spreading, then one year later
plan = pairedsampleplan(grid .|> centroid, control, [-0.003, 0.0, 1.0])

## sample plan figure, deployment geometry

fig = Mke.Figure(size=(500, 400))
ax = Mke.Axis(
    fig[1, 1],
    aspect=1,
    xlabel="Easting [m]",
    ylabel="Northing [m]",
    title="Sample Plan",
)
viz!(ax, treatment, showfacets=true, label="treatment area")
viz!(ax, control, color="gray", showfacets=true, label="control area")
Mke.Legend(fig[1, 2], ax, framevisible=false)
Mke.save(projectdir("demo", "figures", "sample_plan.png"), fig)
fig

##

# 5 cores per composite
samp = CoreSet(plan, 5)

##

# allocate space for the simulation
sim = Simulation((:Ca, :Mg), samp)
# constant feedstock density across all simulations
sim.ρf .= 1e3
# spreading patterns with streaks, as in the tutorial
Q_μ = 3.5
Q = GaussianSimulator(samp, Q_μ)
Q_cov = GaussianCovariance(
    MetricBall((5.0, 500.0)),
    nugget=(3 / 10)^2,
    sill=(3 / 6)^2,
)
# cross-correlated baseline soil concentrations
cs = GaussianCosimulator(samp, (Ca=0.002, Mg=0.001), 0.75)
cs_cov =
    SphericalCovariance(nugget=(0.001 / 10)^2, sill=(0.001 / 6)^2, range=20.0)
# soil spatial structure
ρs = GaussianSimulator(samp, 1000.0)
ρs_cov = SphericalCovariance(nugget=20^2, sill=100^2, range=30.0)
# circular compositing pattern with a radius of 1
stencil = CircleStencil(ncore(sim), 2.0)
# sampler location error with standard deviation of 0.75 m
samplerjitter = Jitter(0.75)
# core location error with standard deviation of 0.02 m (quite precise)
corejitter = Jitter(0.1)
# intentional jitter inside each cell
plannedjitter = GridCentroidJitter(10.0, 5.0, seed=1)
# leaching models for each element
Ca_leaching = ExponentialLeaching(λ=0.4)
Mg_leaching = ExponentialLeaching(λ=0.8)

##

nsim = 10_000

sims = simulationstack(nsim, sim, samp, plan) do
    executeplan!(
        samp,
        plan,
        stencil=stencil,
        plannedjitter=plannedjitter,
        samplerjitter=samplerjitter,
        corejitter=corejitter,
    )

    updategaussian!(Q, samp.points, Q_cov)
    spreading!(rng, sim, plan, Q)

    unmixed!(rng, sim, depth=TriangularDist(0.05, 0.15))

    updategaussian!(ρs, samp.points, ρs_cov)
    rand!(rng, ρs, view(sim.ρs, :))
    #rand!(rng, Normal(1e3, 100), sim.ρs)

    feedstockconcentration!(rng, sim, (Ca=0.07, Mg=0.05), 0.03)

    updategaussian!(cs, samp.points, cs_cov)
    soilconcentration!(rng, sim, cs)

    leaching!(sim, :Ca, Ca_leaching, plan)
    leaching!(sim, :Mg, Mg_leaching, plan)

    massloss!(sim, plan, x -> (x[:Ca] + x[:Mg]) / 2)

    analyze!(rng, sim, (Ca=0.03, Mg=0.03), 0.005)
end

##

tonetcdf(projectdir("demo", "results", "simulations", "demo.nc"), sims)

foreach(
    i -> tocsv(
        projectdir("demo", "results", "simulations", "csv", "demo_$i.csv"),
        sims,
        i,
    ),
    1:nsim,
)

##

fig = Mke.Figure(size=(800, 300))
densegrid = CartesianGrid(minimum(grid), maximum(grid), dims=(150, 150))
γ = Q_cov.γ
r = rand(rng, γ |> GaussianProcess, densegrid, :z => Float64, 6)
for j ∈ 1:3
    i = 1
    ax = Mke.Axis(
        fig[i, j],
        title="Realized Spreading Pattern $(3*(i-1)+j)",
        xlabel="Easting [m]",
        ylabel="Northing [m]",
        aspect=1,
    )
    z = clamp.(r[i+j].z .+ Q_μ, 0.0, Inf)
    viz!(ax, r[i+j].geometry, color=z)
    viz!(ax, control, color=:gray)
end
Mke.Colorbar(
    fig[:, 4],
    limits=(
        clamp(minimum(member -> minimum(member.z) + Q_μ, r), 0.0, Inf),
        maximum(member -> maximum(member.z) + Q_μ, r),
    ),
    colormap=:viridis,
    label="Applied Feedstock [kg/m^2]",
)
Mke.save(projectdir("demo", "figures", "spreading_realizations.png"), fig)
fig

##

fig = Mke.Figure(size=(600, 800))
x = map(p -> p.coords[1], samp.points |> vec)
y = map(p -> p.coords[2], samp.points |> vec)
size = 5

ax = Mke.Axis(
    fig[1, 1],
    aspect=1,
    xlabel="Easting [m]",
    ylabel="Northing [m]",
    title="Application Rate [kg/m^2]",
)
viz!(ax, grid, color=:white, showfacets=true)
s = Mke.scatter!(ax, x, y, color=sim.Q |> vec, markersize=size)
Mke.Colorbar(fig[1, 2], s)

ax = Mke.Axis(
    fig[1, 3],
    aspect=1,
    xlabel="Easting [m]",
    ylabel="Northing [m]",
    title="Soil Bulk Density [kg/m^3]",
)
viz!(ax, grid, color=:white, showfacets=true)
s = Mke.scatter!(
    ax,
    x,
    y,
    color=sim.ρs |> vec,
    colormap=:bamako,
    markersize=size,
)
Mke.Colorbar(fig[1, 4], s)

ax = Mke.Axis(
    fig[2, 1],
    aspect=1,
    xlabel="Easting [m]",
    ylabel="Northing [m]",
    title="Sampling Round (1 - 3)",
)
viz!(ax, grid, color=:white, showfacets=true)
Mke.scatter!(
    ax,
    x,
    y,
    color=repeat(plan.round, 1, 5) |> transpose |> vec,
    colormap=:tab10,
    colorrange=(1, 10),
    markersize=size,
)

ax = Mke.Axis(
    fig[2, 3],
    aspect=1,
    xlabel="Easting [m]",
    ylabel="Northing [m]",
    title="Core Depth [cm]",
)
viz!(ax, grid, color=:white, showfacets=true)
s = Mke.scatter!(
    ax,
    x,
    y,
    color=100 * sim.d |> vec,
    colormap=:reds,
    markersize=size,
)
Mke.Colorbar(fig[2, 4], s)

ax = Mke.Axis(
    fig[3, 1],
    aspect=1,
    xlabel="Easting [m]",
    ylabel="Northing [m]",
    title="Soil Ca Concentration [ppm]",
)
viz!(ax, grid, color=:white, showfacets=true)
s = Mke.scatter!(
    ax,
    x,
    y,
    color=1e6 * sim.cs[:Ca] |> vec,
    colormap=:turbo,
    markersize=size,
)
Mke.Colorbar(fig[3, 2], s)

ax = Mke.Axis(
    fig[3, 3],
    aspect=1,
    xlabel="Easting [m]",
    ylabel="Northing [m]",
    title="Soil Mg Concentration [ppm]",
)
viz!(ax, grid, color=:white, showfacets=true)
s = Mke.scatter!(
    ax,
    x,
    y,
    color=1e6 * sim.cs[:Mg] |> vec,
    colormap=:turbo,
    markersize=size,
)
Mke.Colorbar(fig[3, 4], s)

Mke.save(projectdir("demo", "figures", "simulation.png"), fig)

##
