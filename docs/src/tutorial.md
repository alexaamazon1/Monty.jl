# Simulation Tutorial

Here is an example run through the simulation process for a simple case, from beginning to end.

## Setup

First, we set up the example by loading packages and defining a random number generator for reproducibility. Monty relies heavily on components of [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) for geometric primitives/operations and spatial covariation models.

```@example tutorial
using Monty, Random, Distributions, Meshes
using GeoStatsFunctions, GeoStatsProcesses

import CairoMakie as Mke

const rng = Xoshiro(3)
nothing # hide
```

## Deployment Geometry

Now we define the geometry for our treatment and control areas. They can be anything, but shouldn't overlap. Here I use a simple square, 100 m by 100 m, with a portion of the upper right corner alloted to control samples.
```@example tutorial
# polygon for a treatment area
treatment = PolyArea([(0, 0), (0, 100), (80, 100), (80,50), (100,50), (100,0)])
control = Box(Point(80,50), Point(100,100))

fig = Mke.Figure() # hide
ax = Mke.Axis(fig[1,1], title="Deployment Geometry", xlabel="Easing [m]", ylabel="Northing [m]", aspect=1) # hide
viz!(ax, treatment, label="treatment") # hide
viz!(ax, control, color=:gray, label="control") # hide
fig # hide
```

## Sample Plan

Next, we need to define a sample plan.
* In this case, we'll collect samples at three different times for the imaginary deployment.
* We will also collect the samples at the same locations each time: a *paired* sample plan.
* Even though they are revisited for each *round* of sampling, the sample locations will be randomized.
* Arbitrarily, we'll get 20 samples in the treatment zone and 5 samples in the control zone.
Generating this kind of plan is easy.
```@example tutorial
plan = pairedsampleplan(rng, treatment, 20, control, 5, [-0.03, 0.05, 1.0])
```
The readout indicates that we have 75 total samples over three rounds of sampling: (20 x 3) + (5 x 3). The sample locations are shown below.
```@example tutorial
fig = Mke.Figure() # hide
ax = Mke.Axis(fig[1,1], title="Sample Locations", xlabel="Easing [m]", ylabel="Northing [m]", aspect=1) # hide
viz!(ax, treatment, label="treatment") # hide
viz!(ax, control, color=:gray, label="control") # hide
viz!(ax, plan.points |> PointSet, color=:black, pointsize=12) # hide
fig # hide
```

We can convert the plan into a geo-referenced table using the [`GeoTable`](@ref) constructor, which is reexported by Monty. This is useful for inspection and for exporting the information in a variety of formats (see [`GeoIO.jl`](https://github.com/JuliaEarth/GeoIO.jl)).

```@example tutorial
GeoTable(plan)
```

## Realized Core Locations

Now we execute the sample plan. This means that the target locations are visited, with some error in their locations, and the individual cores for each sample are assigned locations. This is always done using [`executeplan!`](@ref) or [`executeplan`](@ref).

In this case we use 6 cores for each sample, arranged in a "hub and spoke" pattern, meaning a circle with a central point. The core stencil has a radius of 2 meters. We also define some spatial error (jitter) for the sample locations and the individual cores.
```@example tutorial
samp = executeplan(
    plan,
    stencil=HubSpokeStencil(6, 2.0),
    samplerjitter=Jitter(2.0, seed=12),
    corejitter=Jitter(0.1, seed=11)
)
samp # hide
```
Executed plans are stored in a [`CoreSet`](@ref), which is just a matrix of the locations. The readout indicates we now have 450 cores, 6 for each of our 75 samples. Here are the cores, this time with the different colors representing the different sampling rounds.
```@example tutorial
fig = Mke.Figure() # hide
ax = Mke.Axis(fig[1,1], title="Core Locations", xlabel="Easing [m]", ylabel="Northing [m]", aspect=1) # hide
viz!(ax, treatment, label="treatment") # hide
viz!(ax, control, color=:gray, label="control") # hide
viz!(ax, samp.points |> vec |> PointSet, color=repeat(plan.round, inner=ncore(samp)), pointsize=6) # hide
fig # hide
```

The realized core locations are not associated with the rest of the metadata in the sample plan, but we can make a [`GeoTable`](@ref) by calling the constructor with both groups of information. This time the table includes the "core" column, which indicates the core number for each sample location.

```@example tutorial
GeoTable(samp, plan)
```

## Simulation

We have all the information about sample locations and times that we need, in addition to whether samples are controls.

To simulate, we define a container that makes it easy to keep track of simulated quantities and set their values. This is the [`Simulation`](@ref) type. Because it's just a container, we only need the number of samples and cores per sample to allocate one.

Below a simulation is started with `Ca`, `Mg`, and `Zr` as the target analytes/elements. By passing `samp` to the constructor, the simulation is created with the correct size automatically.

```@example tutorial
sim = Simulation((:Ca, :Mg, :Zr), samp)
```
The readout for a [`Simulation`](@ref) is long, but it's an exhaustive summary of the information inside the container and it indicates if fields have not been defined yet. Undefined fields have `NaN` values in their summaries. Initially, everything except the core cross-sectional area is undefined. This area is usually not important so it's assigned a default value, but it can also be declared.

Now we have to fill in the parameters of a simulation,  defining each of them for every core in the deployment.

#### Application Rate

First, applied feedstock. I'll assume an application of, on average, 3 kg/m$^2$. I ignore moisture for now. Feedstock is generally applied by spreaders hitched to tractors, which are driven over fields in linear passes. As such, we might assume the variability in application rate over our treatment area has some spatial structure. Maybe it has a vaguely striped and streaky pattern.

This can be represented using an appropriate variogram, which defines correlations between points based on their distances from each other. A lot can be said about variograms, but here I will just demonstrate how to use one and use it to fill in the application rate for the simulation. Instead of a [`Variogram`](https://juliaearth.github.io/GeoStatsDocs/stable/variography/theoretical.html) I define a `Covariance` which is more convenient and what the Monty functions expect.

```@example tutorial
# average application rate
Q = 3.0
# this "ball" defines the range on each axis
ball = MetricBall((5.0, 1000.0))
# spatial structure
C = CubicCovariance(
    ball,
    nugget=(0.1Q)^2, # 10 % nugget effect
    sill=(0.25Q)^2,  # 25 % overall standard deviation
)
```
What does this structure look like? It's helpful to realize the spatial structure on a dense grid, to get a sense for it. One realization is shown below. The overall structure is streaky with some noise over very short distances from the nugget effect.
```@example tutorial
gp = CubicVariogram(ball, nugget=(0.1Q)^2, sill=(0.25Q)^2) |> GaussianProcess # hide
grid = CartesianGrid((0.0, 0.0), (100.0, 100.0), dims=(200, 200)) # hide
r = rand(rng, gp, grid, :z => Float64) # hide
fig = Mke.Figure() # hide
ax = Mke.Axis(fig[1,1], title="Example Spreading Pattern", xlabel="Easing [m]", ylabel="Northing [m]", aspect=1) # hide
viz!(ax, r.geometry, color=r.z) # hide
fig # hide
```

From the covariance type defined above, Monty will assemble a multivariate Gaussian distribution with the proper covariance matrix. This distribution can then be used to fill in the application rate for the simulation once or, for repeated simulations, as many times as desired.

The constructor needs our [`CoreSet`](@ref), which we named `samp`, because `samp` contains all the realized core locations. There is also a [`spreading!`](@ref) function that draws from the distribution and puts the values into the correct places, *ignoring control points automatically*.
```@example tutorial
X = Monty.MvNormal(Q, C, samp)
spreading!(rng, sim, plan, X)
sim
```
Notice that the readout for the simulation has been updated. It shows a non-`NaN` summary for `Q`, which is our symbol for the application rate. The summary shows the mean value, the standard deviation, and the extrema. Note that the mean value is not 3, as expected. This is because of the control samples. Of our 25 sample locations, 20 are treatment and 5 are control. The expected average is 2.4, close to what we see above.

#### Pre-Spreading Soil Concentrations

It's reasonable to expect the baseline soil, before any feedstock is spread, to have some spatial structure as well. We expect calcium, magnesium, and zirconium concentrations in the soil to be similar to each other over short distances. Further, we might expect calcium and magnesium to be cross-correlated, meaning we tend to find high calcium wherever we find high magnesium, and vice versa.

This can also be represented in Monty. For two cross-correlated fields like `Ca` and `Mg`, Monty has a [`GaussianCosimulator`](@ref) type that handles spatial autocorrelation in two analytes and cross-correlation between them.
```@example tutorial
# mean values for Ca and Mg in soil
μ = (Ca=0.002, Mg=0.001)
c = SphericalCovariance(
    nugget=1e-4^2,
    sill=3e-4^2,
    range=20.0
)
ρ = 0.7
gc = GaussianCosimulator(samp, μ, c, ρ)
```
This type also has a convenience function associated with it, to put the realization into the simulation.
```@example tutorial
soilconcentration!(rng, sim, gc)
```
Here is what the baseline soil concentrations actually look like:
```@example tutorial
fig = Mke.Figure(size=(800,400)) # hide
ax = Mke.Axis(fig[1, 1], aspect=1, title="Ca", xlabel="Easing [m]", ylabel="Northing [m]") # hide
viz!(ax, treatment, alpha=0.3) # hide
viz!(ax, control, color=:gray, alpha=0.3) # hide
viz!(ax, samp.points |> vec |> PointSet, color=sim.cs[:Ca] |> vec, pointsize=8) # hide
ax = Mke.Axis(fig[1, 2], aspect=1, title="Mg", xlabel="Easing [m]", ylabel="Northing [m]") # hide
viz!(ax, treatment, alpha=0.3) # hide
viz!(ax, control, color=:gray, alpha=0.3) # hide
viz!(ax, samp.points |> vec |> PointSet, color=sim.cs[:Mg] |> vec, pointsize=8) # hide
fig # hide
```
Finally, we assume no spatial structure for `Zr`, just to keep things moving.
```@example tutorial
μ = 10e-6 # 10 ppm on average
X = Normal(μ, μ/5)
soilconcentration!(rng, sim, :Zr, X)
```
The simulation readout is getting filled in.
```@example tutorial
sim # hide
```

#### Sample Depth & Mixing

The sample depth and the mixing profile of feedstock in soil interact. In general, deeper samples contain more soil per unit feedstock because most of the feedstock is near the surface. But it's also possible for shallow samples to collect only part of the original applied feedstock because it has been mixed more deeply than the sample depth.

Monty handles sample depth and mixing profile jointly. The sample depth can be defined by a non-negative probability distribution and the mixing profile by a distribution over depth *for each core*. For example, imagine the feedstock mixing profile looks like a wedge, tapering off to zero at some depth which varies naturally from core to core. The sample depth also varies independently at each location.
```@example tutorial
triangularmixing!(
    rng,
    sim,
    depth=truncated(Normal(0.125, 0.025), 0.0, Inf),
    upper=Exponential(0.05)
)

```
This operation sets the depth (`d`) and the feedstock fraction (`γ`) for all cores in the simulation. Here is the updated read-out.
```@example tutorial
sim
```

#### Leaching & Mass Loss

Just a couple more quantities worth attention before generating a synthetic dataset.

The amount of an element left in feedstock as time proceeds is defined by a "leaching model" in Monty. These models simply define functions that are
* equal to one for time < 0
* equal to zero for time = 0
* monotonically increasing after time = 0

For example, here is a leaching model that represents exponential loss of an element over time, with a seasonal fluctuation on top. The value of the line represents the *fraction of feedstock present/lost*. Before time zero, this is always equal to one because no feedstock has been spread yet.
```@example tutorial
w = SeasonalLeaching(λ=2.0, C=0.9, power=1, phase=π / 2)
t = LinRange(-1.0, 3.0, 1000)
fig = Mke.Figure()
ax = Mke.Axis(
    fig[1, 1],
    title="Seasonal Leaching Model",
    xlabel="Time [years]",
    ylabel="Feedstock Loss Fraction [-]"
)
Mke.lines!(ax, t, w.(t))
fig # hide
```
We can apply these models to mobile cations in the simulation. Immobile elements must also be assigned leaching models (no leaching).
```@example tutorial
leaching!(sim, :Ca, ExponentialLeaching(λ=0.8), plan)
leaching!(sim, :Mg, ExponentialLeaching(λ=1.2), plan)
leaching!(sim, :Zr, NoLeaching(), plan)
```

Additionally, we can define the mass loss fraction for all the feedstock, not just elements that are leaching out of it. For example, the material may have lost 30% of a mobile element, but 50 % of its total mass because it's composed of other elements like silicon and oxygen that are also liberated.

The total mass loss can be any function of the individual elemental loss fractions that produces a fractional value between zero and one. It should make physical sense, but Monty will not check that for you.

To keep it simple, I'll just use an average of the Ca and Mg loss fractions.
```@example tutorial
massloss!(sim, plan, x -> (x[:Ca] / 2 + x[:Mg] / 2))
sim # hide
```

#### Other Physical Quantities

The rest of the parameters are uncomplicated. I'll assign them independently with appropriate distributions.
```@example tutorial
# feedstock bulk density
sim.ρf .= 2e3
# soil bulk density
rand!(rng, Normal(1e3, 100), sim.ρs)
# feedstock concentrations with 2 % variability
feedstockconcentration!(rng, sim, (Ca=0.07, Mg=0.06, Zr=0.001), 0.02)
```
Now the simulation has complete information. Here's the readout.
```@example tutorial
sim # hide
```

## Generating Measurements

Once all the information is in place, the rest is very easy. There are three simple steps
1. apply the mixing model to each core, using the values defined in the simulation
2. composite the cores for each sample (by summing them)
3. "measure" the composite samples by applying analytical noise to the concentrations and the sample masses

The simulation already has memory allocated for these operations, so they happen internally.
```@example tutorial
core!(sim)
sim.cores[1,1]
```
```@example tutorial
composite!(sim)
sim.composites[1]
```
Here I use measurement noise of 4 %, 5 %, and 7 %, representing error through the whole process of splitting, preparing, and analyzing the samples. I also use 0.5 % error for the mass.
```@example tutorial
measure!(rng, sim, (Ca=0.04, Mg=0.05, Zr=0.07), 0.005)
sim.measurements[1]
```
Alternatively, all of this can be executed at once with the [`analyze!`](@ref) function, which calls [`core!`](@ref), [`composite!`](@ref), and [`measure!`](@ref) for you.
```@example tutorial
analyze!(rng, sim, (Ca=0.04, Mg=0.05, Zr=0.07), 0.005)
sim.measurements[1]
```

The final, simulated measurements are all in the `sim.measurements` vector. They can be arranged into a [`GeoTable`](@ref), and then exported in whatever format is convenient (again, see [`GeoIO.jl`](https://github.com/JuliaEarth/GeoIO.jl)).
```@example tutorial
table = GeoTable(sim, samp, plan)
```
The `geometry` column above represents an average of the realized core locations. To use only the planned locations, omit `samp` from the constructor with `GeoTable(sim, plan)`.

## Analysis

With a synthetic dataset in hand, we can do whatever we want with it. In many cases, much of the process above will be repeated to generate a large number of datasets using the same variability and physical properties. After defining a sample plan, we can realize lots of core sets and then generate lots of hypothetical samples and measurements. There are a few more tricks to run repeated simulations efficiently, with zero memory overhead, using the in-place methods highlighted above. An example is also shown on the [Efficient Simulation](efficiency.md) page.

Here, however, I'll just do a basic visualization of the individual dataset generated above. First, it can be converted to a `DataFrame`.
```@example tutorial
using DataFrames

df = table |> DataFrame
df[1:5,:]
```
Here is a summary of the concentrations in each round, converted to units of ppm by mass.
```@example tutorial
fig = Mke.Figure(size=(1000,400)) # hide
for (i,a) in enumerate([:Ca, :Mg, :Zr]) # hide
    ax = Mke.Axis(fig[1, i], title="$a", xlabel="Sampling Round", ylabel="Simulated $a Concentration [ppm]") # hide
    noise = randn(size(df,1)) / 20 # hide
    control = df[!,:control] .|> Int # hide
    dodge = (control .- 0.5) ./ 3 # hide
    Mke.scatter!(ax, df[!,:round] .+ noise .+ dodge, 1e6 * df[!,a], color=control .+ 1, strokewidth=0.25, colormap=:tab10, colorrange=(1,10)) # hide
end # hide
fig # hide
```
The blue dots show the 20 treatment points for each round and the orange ones show the 5 control points. The horizontal axis is the sampling round. This is easier to visualize than the actual time of each sampling round, because the first and second rounds occur almost at the same time (before and after spreading). A few cursory observations:
* For each element, there is no clear difference in the concentrations for control points over the different sampling rounds. All we see is the natural variation of the baseline soil, averaged into our composite core for each sample and then measured with analytical error.
* Concentrations for Ca and Mg go up after spreading, then down again. This is the expected pattern for spreading then weathering and leaching.
* Concentrations for Zr, however, go up and stay up. This is because we used a [`NoLeaching`](@ref) model for Zr. We assumed that it was perfectly immobile.
* Although there are relatively few control points, the *spread* in the points seems bigger for the treatment points. This is a consequence of the variability in spreading/application of feedstock, which we designed into the simulation above.
* This is an example of nice, well-behaved data. It would be considered high-quality MRV data for a real deployment. In general, real data will be noisier and messier.
