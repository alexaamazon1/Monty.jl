# Big Simulation & Analysis

Here is an example of a single simulation with a really large number of samples over four different rounds of sampling. The example shows how to use 32-bit floats for the whole simulation, saving some memory, with some basic plotting at the end. The explanation is more brief here compared to the [tutorial](tutorial.md) and other examples.

## Setup

Load the necessary packages and seed a random number generator.

```@example bigsim
using Monty
using Random
using Distributions
using GeoStatsFunctions
using DataFrames
using Meshes

import CairoMakie as Mke

const rng = Xoshiro()
nothing # hide
```

## Sampling Plan

Here, the field geometry is just a big circle, roughly equivalent to 1000 acres. We take 500 samples in 4 different sampling rounds, each sample with 6 compsited soil cores. In each round, the samples are randomly located. We only simulate calcium and magnesium. We also neglect control samples for this example.

```@example bigsim
field = Ball((0f0, 0f0), 1135f0)
times = Float32[0, 1, 2, 3]
plan = randomsampleplan(rng, field, 500, times)
samp = CoreSet(plan, 6)
sim = Simulation((:Ca, :Mg), samp)
fig = Mke.Figure(size=(500, 500)) # hide
ax = Mke.Axis( # hide
    fig[1, 1], # hide
    xlabel="Easting [m]", # hide
    ylabel="Northing [m]", # hide
    title="Sample Locations", # hide
    aspect=1, # hide
) # hide
viz!(ax, field, color="gray", alpha=0.4) # hide
Mke.scatter!( # hide
    ax, # hide
    map(p -> p.coords[1], plan.points), # hide
    map(p -> p.coords[2], plan.points), # hide
    color=plan.round, # hide
    markersize=5, # hide
) # hide
fig # hide
```

The plot above gives a sense for how dense the sampling is. The colors represent the 4 different sampling rounds.

## Define Model Parameters

Here we initialize types for representing spatially correlated fields in the simulation for the application rate and the baseline soil concentration. We also define the sampling jitter, compositing stencil, and leaching models (for both elements).

```@example bigsim
Q = GaussianSimulator(samp, 4.0f0)
Q_cov = GaussianCovariance(
    MetricBall((2.0f1, 5.0f2)),
    nugget=(4.0f-1)^2,
    sill=8.0f-1,
)

cs = GaussianCosimulator(samp, (Ca=3.0f-3, Mg=1.0f-3), 5.0f-1)
cs_cov = SphericalCovariance(
    nugget=Float32(5e-5)^2,
    sill=Float32(1e-4)^2,
    range=20.0f0,
)

stencil = CircleStencil(ncore(sim), 3.0f0)
samplerjitter = Jitter(5.0f0)
corejitter = Jitter(1.0f-1)

Ca_leaching = ExponentialLeaching(λ=0.4f0, C=1.0f0, σ=0.0f0)
Mg_leaching = ExponentialLeaching(λ=0.6f0, C=1.0f0, σ=0.0f0)
nothing # hide
```

## Simulate

Now we execute the sampling plan and fill in all the fields of our simulation. Note (in the print out after the code block) that the `Simulation` has element type `Float32`. With a very large number of cores, the covariance matrix for cross-correlated fields is very large. In this case the `GaussianCosimulator` contains a 24,000 x 24,000 matrix for two fields (Ca and Mg) over 12,000 individual cores. It might be fine to use 64-bit floats, but we use 32-bit just as an example of how to save some memory, in case you have a small computer or an even larger number of cores.

```@example bigsim
executeplan!(
    samp,
    plan,
    stencil=stencil,
    samplerjitter=samplerjitter,
    corejitter=corejitter,
)

updategaussian!(Q, samp.points, Q_cov)
spreading!(rng, sim, plan, Q)

updategaussian!(cs, samp.points, cs_cov)
soilconcentration!(rng, sim, cs)

sim.ρf .= 1.0f3

exponentialmixing!(
    rng,
    sim,
    depth=TriangularDist(0.08f0, 0.12f0),
    scale=Uniform(1.0f-2, 3.0f-2),
)

rand!(rng, Normal(1.0f3, 5.0f1), sim.ρs)

feedstockconcentration!(rng, sim, (Ca=0.07f0, Mg=0.05f0), 0.05f0)

rand!(rng, Normal(0.002, 0.002 / 10), sim.cs[:Ca])
rand!(rng, Normal(0.001, 0.001 / 10), sim.cs[:Mg])

leaching!(sim, :Ca, Ca_leaching, plan)
leaching!(sim, :Mg, Mg_leaching, plan)

massloss!(sim, plan, x -> (x[:Ca] / 2 + x[:Mg] / 2))

analyze!(rng, sim, (Ca=0.05f0, Mg=0.06f0), 0.005f0)
sim #hide
```

# Results

First, we convert the simulation to a `GeoTable` and then a `DataFrame`, for convenience. Here are the first five rows.
```@example bigsim
df = GeoTable(sim, plan) |> DataFrame
df[1:5,:] # hide
```

Then we look at box plots of each sampling round's concentration. Because we have so many samples, the box plots are beautiful and and the mean values fall nicely along an exponential curve. Also, note that the variability of both concentrations is smaller as time goes on. This is a consequence of the variability in application rates. As feedstock is dissolved and leached out of the mixture, this variability is attenuated.

```@example bigsim
fig = Mke.Figure(size=(800, 400)) # hide
ax1 = Mke.Axis( # hide
    fig[1, 1], # hide
    title="Calcium", # hide
    ylabel="Concentration [ppm]", # hide
    xlabel="Time [yr]", # hide
) # hide
Mke.boxplot!(ax1, df.time, 1e6 * df.Ca) # hide
ax2 = Mke.Axis( # hide
    fig[1, 2], # hide
    title="Magnesium", # hide
    ylabel="Concentration [ppm]", # hide
    xlabel="Time [yr]", # hide
) # hide
Mke.boxplot!(ax2, df.time, 1e6 * df.Mg) # hide
Mke.linkaxes!(ax1, ax2) # hide
fig # hide
```
