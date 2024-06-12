# Efficient Simulation

The [tutorial simulation](tutorial.md) walks through the steps of a single simulation. In many cases, we will want to run a large number of *realizations* for a hypothetical deployment, sampling from all the variable inputs to see how much the results vary. This basic Monty Carly approach is flexible and can be used to answer lots of questions (hence the name `Monty`). But, the quality of Monty Carlo results depends on the number of samples: more is better. Setting up efficient repeated simulations requires some attention and depends almost entirely on whether iterations avoid unnecessary memory allocation. The example below shows how to do this.

## Setup

First, we load the modules

```@example efficient
using Monty, Random, Distributions, Meshes
using DataFrames
using GeoStatsFunctions

import CairoMakie as Mke

const rng = Xoshiro(42)
nothing # hide
```
and set up the geometry for a hypothetical deployment. For this example, I use a simple rectangular grid split up into alternating treatment/control strips. The grid has 100 total cells, each of which will have one composite sample taken from it.
```@example efficient
grid = CartesianGrid((10,10), (0.0, 0.0), (10.0, 10.0))
treatment = reduce(vcat, [grid[i,:] for i ∈ 1:2:10])
control = reduce(vcat, [grid[i,:] for i ∈ 2:2:10])
treatment_points = treatment .|> centroid
control_points = control .|> centroid

fig = Mke.Figure() # hide
ax = Mke.Axis(fig[1,1], xlabel="Easing [m]", ylabel="Northing [m]") # hide
viz!(ax, treatment) # hide
viz!(ax, control, color=:gray) # hide
viz!(ax, grid .|> centroid, color=:black, pointsize=8) # hide
fig # hide
```
In this example, we use an intentional "grid centroid jitter," which means every realization of the sample plan moves the sample points to random locations inside each grid cell. This can help avoid potential issues with putting sample points along straight lines. Here is a quick visualization of one realization of this jitter on the grid.
```@example efficient
plannedjitter = GridCentroidJitter(10.0, 2.0)
points = plannedjitter.(grid .|> centroid)

fig = Mke.Figure() # hide
ax = Mke.Axis(fig[1,1], xlabel="Easing [m]", ylabel="Northing [m]") # hide
viz!(ax, grid, color=:white, facetcolor=:gray, showfacets=true) # hide
viz!(ax, points, color=:black, pointsize=8) # hide
fig # hide
```

## Pre-allocation

Running efficient simulations means avoiding unnecessary memory allocations inside the loop. To avoid them, we have to allocate the memory we need for simulations once and use it for each realization, shuffling the final output of each realization into some kind of data structure along the way and returning those results. We also need to set the parameters and distributions that each simulation is realized from.

Here I allocate everything needed for repeated simulations. First, I create the sample `plan` and allocate space for cores (`samp`) and simulations `sim`. Then I define the variables needed to realize an individual simulation, which captures both operational and physical/chemical variability.
```@example efficient
# location-paired sampling at t=0 and t=1
plan = pairedsampleplan(treatment_points, control_points, [0.0, 1.0])
# allocate space for realized core locations, 10 composites per sample
samp = CoreSet(plan, 10)
# allocate space for the simulation
sim = Simulation((:Ca, :Mg), samp)

# spreading patterns with streaks, as in the tutorial
Q = GaussianSimulator(samp, 4.0)
Q_cov = GaussianCovariance(
    MetricBall((10.0, 1000.0)),
    nugget=(4.0 / 10)^2,
    sill=(4.0 / 3)^2,
)
# circular compositing pattern with a radius of 1
stencil = CircleStencil(ncore(sim), 1.0)
# sampler location error with standard deviation of 2
samplerjitter = Jitter(2.0)
# core location error with standard deviation of 0.1
corejitter = Jitter(0.1)

# constant feedstock density across all simulations
sim.ρf .= 3e3
# cross-correlated baseline soil concentrations
cs = GaussianCosimulator(samp, (Ca=0.002, Mg=0.001), 0.7)
cs_cov = SphericalCovariance(nugget=5e-5^2, sill=1e-4^2, range=20.0)
# variable sample depths
depth = TriangularDist(0.06, 0.14)
# somewhat variable soil bulk density
ρs = Normal(1e3, 50)
# feedstock mean concentrations with 5 % variability
cf = (Ca=0.07, Mg=0.05)
σf = 0.05
# leaching models for each element
Ca_leaching = ExponentialLeaching(λ=1.0)
Mg_leaching = ExponentialLeaching(λ=1.5)
# analytical noise
σ_analysis = (Ca=0.04, Mg=0.03)
nothing # hide
```

## Simulation Stacks

Monty has a convenience function for repeating a simulation and dumping the results into a final stack of arrays. This function is meant to be used with Julia's [do-block syntax](https://docs.julialang.org/en/v1/manual/functions/#Do-Block-Syntax-for-Function-Arguments). This means that we write the simulation routine into the function call directly. See below.
```@example efficient
sims = simulationstack(100, sim, samp, plan) do

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
    clamp!(sim.Q, 0.0, Inf)

    unmixed!(rng, sim, depth=depth)

    rand!(rng, ρs, sim.ρs)

    feedstockconcentration!(rng, sim, cf, σf)

    updategaussian!(cs, samp.points, cs_cov)
    soilconcentration!(rng, sim, cs)
    clamp!(sim.cs[:Ca], 0.0, Inf)
    clamp!(sim.cs[:Mg], 0.0, Inf)

    leaching!(sim, :Ca, Ca_leaching, plan)
    leaching!(sim, :Mg, Mg_leaching, plan)

    massloss!(sim, plan, x -> (x[:Ca] / 2 + x[:Mg] / 2))

    analyze!(rng, sim, (Ca=0.01, Mg=0.005), 0.005)

end
sims # hide
```
The `simulationstack` function above is passed `100`, the number of realizations, then the data structures we need to simulate (`sim`, `samp`, `plan`). Then the `do` block contains the entire script for our simulation.

The `simulationstack` function handles the storage of results automatically, returning a stack of arrays in the form of a `DimStack` from [`DimensionalData.jl`](https://github.com/rafaqz/DimensionalData.jl). As the read-out above indicates, the simulation results are packed into the `data` member of the results, with associated information like control flags in the other members.

For example, to see the first realization, we can just slice it out of the results:
```@example efficient
sims[:data][1,:,:]
```
What happens next depends on what we want to know. If you just want to save the data, use [`tonetcdf`](@ref) or [`tocsv`](@ref). You can also convert the entire body of simulations directly into a `DataFrame`. Here are the first 5 rows of the dataframe for all simulations:
```@example efficient
df = sims |> DataFrame
df[1:5,:] # hide
```
Or you can slice/select from the simulations and convert *that* portion into a frame. Here are the first 5 rows of this dataframe for only the first realization:
```@example efficient
df = unstack(sims[1,:,:] |> DataFrame, :analyte, :data)
df[1:5,:] # hide
```

## Notes

* In the example above, *almost all* of the computation time is spent updating the Gaussian simulators (performing [`cholesky!`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.cholesky!) factorizations). Without spatial correlation/structure, the simulations are very fast. Benchmarks with a similar simulation script indicate that simulations without spatial correlation are about 1000 times faster. This is only a very rough comparison though.
    * Without spatial correlation, the computational cost scales linearly. These simulations will be extremely fast, on the order of tens of *microseconds*. So, it would take about 1 second to run 100k realizations of a hypothetical deployment without spatial correlation.
    * With spatial correlation, the timing depends roughly *cubically* on the number of total cores in the simulation because of the [Cholesky factorizations](https://en.wikipedia.org/wiki/Cholesky_decomposition). So, it may vary a lot. To speed things up, one option is to use spatial correlation but not relocate the sample points for every realization. This avoids refactorizing the covariance matrix and the values of the points can still be drawn random from a persistent multivariate normal distribution.
* Benchmarks also verify that the simulation loop above allocates zero memory, which is good.
* What goes into the simulation script, inside the `do` block, is up to you. If any parts of the simulation are shared across realizations, leave them out of the script. For example, above the `ρf` field of the simulation is assumed to be the same in all realizations, so it's set outside the simulation loop.
* I think the `do`-block syntax helps keep the focus on simulation choices, instead of how the results arrays are allocated and filled. We just have to write the script, and then the looping and storing for each realization is handled out of view.
