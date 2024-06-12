# Reference

## Sample Planning

### Generating Points

```@docs
randin
randin!
gridpointoverlay
jittergrid
```

### Sample Plans

```@docs
SamplePlan
nsample
boundingbox
pairedsampleplan
randomsampleplan
```

## Location Jitter

Monty has three "jitter" types meant to be used as callable objects or "functors." This means that a jitter object is created with some parameters, then used like a function and applied to sample/core locations. Each jitter type can be called on a `Point` and it will return a jittered `Point`.

```@docs
jitterpoint
jitterpoints!
NoJitter
Jitter
GridCentroidJitter
```

## Core Stencils

Like the jitter types above, core stencil types are all used as callable objects or "functors." This means that a stencil is created with some parameters, then used like a function to generate core locations. Generally, stencils don't have to be used directly in Monty. They are passed to [`executeplan!`](@ref), which handles the function calling. To use one directly, however, just requires passing a location and core number to the object

    stencil(point::Point, n::Integer)

which returns the core location as another `Point`.

```@docs
SingleCoreStencil
RandomStencil
CircleStencil
HubSpokeStencil
LineStencil
```

## Plan Execution

```@docs
CoreSet
executeplan
executeplan!
```

## Mixing Model

```@docs
feedstockfraction
feedstockmass
feedstockelementmass
feedstockthickness
soilmass
soilelementmass
samplemass
Sample
nanalyte
analytes
mixing
```

## Weathering & Leaching Models

Leaching models define the fraction of a feedstock's element *lost* at any point in time. Each leaching model is a scalar function of time. They always return values that are
* fractional, meaning $\in [0,1]$
* exactly 1 for negative time
* exactly zero at time zero
* monotonically increasing after time zero (no reverse leaching)

```@docs
NoLeaching
ExponentialLeaching
MultiExponentialLeaching
SeasonalLeaching
```

Here are some a example leaching models:
```@example
using Monty # hide
import CairoMakie as Mke # hide
# time points where models are evaluated
t = LinRange(-1.0, 4.0, 51) |> collect
# models
exp_model = ExponentialLeaching(λ=1.5, C=0.9)
multi_exp_model = MultiExponentialLeaching(λ=(10.0, 1.0), C=0.9)
seasonal_model = SeasonalLeaching(λ=1.0, floor=0.05, power=2, C=0.85)

fig = Mke.Figure(size=(800,400)) # hide
ax = Mke.Axis(fig[1,1], xlabel="time [yr]", ylabel="leached fraction") # hide
Mke.scatter!(ax, t, exp_model.(t), alpha=0.7, label="exp_model") # hide
Mke.scatter!(ax, t, multi_exp_model.(t), alpha=0.7, label="multi_exp_model") # hide
Mke.scatter!(ax, t, seasonal_model.(t), alpha=0.7, label="seasonal_model") # hide
Mke.Legend(fig[1,2], ax, framevisible=false) # hide
fig # hide
```

## Simulation

The main simulation function is [`simulationstack`](@ref) and an example is shown in the [Efficient Simulations](efficiency.md) page.

```@docs
moisturefraction
moistureratio
spreading!
unmixed!
triangularmixing!
uniformmixing!
exponentialmixing!
feedstockconcentration!
soilconcentration!
leaching!
massloss!
core!
composite!
measure!
analyze!
simulationstack
tonetcdf
tocsv
```

### Container

```@docs
Monty.Simulating.Simulation
clear!
```

### Gaussian Simulators

There are two types for simulating from multivariate gaussian distributions, which are useful for generating vectors with correlated elements. Both of them can be updated from a `Covariance` model using the `updategaussian!` methods, which recomputes the necessary Cholesky factorization in-place using the covariance model and a group of points. Both types can also be sampled with `rand` and `rand!`.

```@docs
GaussianSimulator
GaussianCosimulator
updategaussian!
Monty.Gaussian.rand!
Monty.Gaussian.rand
```

## Measurement

```@docs
noise
Measurement
measure
```

## Tabulation

Each of the main composite types in Monty can be exported to a [`GeoTable`](https://juliaearth.github.io/GeoStatsDocs/stable/data.html), which can then be exported to a number of file formats using [`GeoIO.jl`](https://github.com/JuliaEarth/GeoIO.jl).

```@docs
GeoTable
```

## CDR Potential

There are a few convenience functions for computing the CDR potential of a feedstock.

```@docs
cationCO2
cdrpotential
```
