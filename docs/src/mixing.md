# Mixing Tutorial

The central element of `Monty` is a mixing model describing how the concentration of an individual soil sample is determined by a mixture of soil and (possibly weathered) feedstock. This component of the package is deterministic. For a single core, physical parameters for the feedstock and the soil are passed in and the concentrations of specified elements are computed, in addition to the core's mass.

The mixing model assumes that baseline soil is vertically homogenous and that the weathering/dissolution of feedstock is also vertically homogenous. It does not, however, assume that the *mixing* of feedstock into the soil is vertically homogenous. It allows for different mixing profiles to be prescribed.

## Parameters

The inputs for the mixing model are tabulated below. These represent the physical and chemical conditions in a single soil core.

| Parameter | SI Unit | Description |
| :-------: | :-----: | :---------- |
| `γ` | kg/kg | fraction of applied feedstock present in the sampling depth |
| `d` | m | sample depth |
| `a` | m$^2$ | cross-sectional area |
| `Q` | kg/m$^2$ | feedstock application rate (dry material) |
| `ρf` | kg/m$^3$ | feedstock bulk density |
| `cf` | kg/kg | elemental concentration(s) in feedstock |
| `ρs` | kg/m$^3$ | soil bulk density |
| `cs` | kg/kg | elemental concentration(s) in soil |
| `𝓁` | kg/kg | loss fraction(s) of elemental concentrations in feedstock due to weathering |
| `ℒ` | kg/kg | loss fraction of total feedstock mass due to weathering |

## Examples

First, load the relevant packages.

```@example mixing
using Unitful
using Monty
using Distributions
using BenchmarkTools
```

To pull a hypothetical soil sample, use the [`mixing`](@ref) function. For this example, we use [unitful](https://painterqubits.github.io/Unitful.jl/stable/) numbers, which automatically check that the dimensions of the inputs and outputs are sensible.

```@example mixing
mixing(
    1.0u"kg/kg",
    0.1u"m",
    1.25e-3u"m^2",
    2.2u"kg/m^2",
    2e3u"kg/m^3",
    5e4u"ppm",
    900.0u"kg/m^3",
    1e3u"ppm",
    0.0u"kg/kg",
    0.0u"kg/kg",
)
```

Calling the mixing model creates a [`Sample`](@ref), which is printed out above with information about its type and its contents. The numeric type of the sample is complex because we used units. The value for the `mass` is straightforward. Because we didn't specify the name of the analyte, the default name is simply `analyte` and its concentration is printed with the same units as the inputs. The hypothetical soil core has a concentration of about 2182 ppm for this arbitrary analyte.

If we take out the units, it's a little simpler. Without units, the concentrations are assumed to be mass fractions (kg/kg) and the `Sample` type is simply `Float64`.
```@example mixing
mixing(1.0, 0.1, 1.25e-3, 2.2, 2e3, 0.05, 900.0, 0.001, 0.0, 0.0)
```

Most of the time, we're interested in multiple elements. The mixing model runs multiple elements automatically if the concentrations are specified by [named tuples](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple). This is a fast and simple way for the mixing model to keep track of different analyte names and it forces us to be specific. For example, to use two elements, just specify their concentrations in both feedstock and soil, as well as their loss fractions (`𝓁`) from feedstock.
```@example mixing
cf = (Ca=0.05, Mg=0.06)
cs = (Ca=0.001, Mg=0.0005)
𝓁 = (Ca=0.0, Mg=0.0)
mixing(1.0, 0.1, 1.25e-3, 2.2, 2e3, cf, 900.0, cs, 𝓁, 0.0)
```

!!! warning
    The model is very strict about using consistent numeric types in the named tuples. You can't use `(Ca=0.0, Mg=0)`, for example, because it mixes an integer with a float.

Instead of supplying the feedstock mixing fraction `γ` directly, we can supply the entire vertical mixing profile `Γ`, which must be a univariate, non-negative probability distribution. Specifically, it must be a [`UnivariateDistribution`](https://juliastats.org/Distributions.jl/stable/univariate/) for which `minimum(Γ) >= 0`.
```@example mixing
Γ = truncated(Normal(0.0, 0.05), 0.0, Inf) # a half normal with σ = 5 cm
cf = (Ca=0.05, Mg=0.06)
cs = (Ca=0.001, Mg=0.0005)
𝓁 = (Ca=0.0, Mg=0.0)
mixing(Γ, 0.1, 1.25e-3, 2.2, 2e3, cf, 900.0, cs, 𝓁, 0.0)
```

Finally, the mixing model is fast. Even though it creates [`Sample`](@ref) instances, these are mainly just wrappers for named tuples. Here is a benchmark.
```@example mixing
cf = (Ca=0.05, Mg=0.06)
cs = (Ca=0.001, Mg=0.0005)
𝓁 = (Ca=0.0, Mg=0.0)
@benchmark mixing(1.0, 0.1, 1.25e-3, 2.2, 2e3, $cf, 900.0, $cs, $𝓁, 0.0)
```
The model is simple, so it ought to be fast, but it's good confirmation that it allocates no memory and takes only tens of nanoseconds to run.

For one last example, we can also use multiple analytes with unitful numbers. No problem.
```@example mixing
mixing(
    1.0u"kg/kg",
    0.1u"m",
    1.25e-3u"m^2",
    2.2u"kg/m^2",
    2e3u"kg/m^3",
    (Ca=5e4u"ppm", Mg=6e4u"ppm"),
    900.0u"kg/m^3",
    (Ca=1e3u"ppm", Mg=1200.0u"ppm"),
    (Ca=0.0u"kg/kg", Mg=0.0u"kg/kg"),
    0.0u"kg/kg",
)
```

## Combining samples

The physical combination of cores and samples (compositing) is easy. Just add them together with `+`, use `sum`, or any other version of addition.
```@example mixing
a = mixing(1.0, 0.11, 1.25e-3, 2.2, 2e3, 0.05, 900.0, 0.001, 0.0, 0.0)
a # hide
```
```@example mixing
b = mixing(1.0, 0.09, 1.25e-3, 2.2, 2e3, 0.05, 900.0, 0.001, 0.0, 0.0)
b # hide
```
```@example mixing
c = a + b
c # hide
```
Adding the two individual cores gives you a [`Sample`](@ref) with 2 cores and you can verify that the concentration is a mass-weighted average of the two individual concentrations. There's no limit to the number of cores you can combine.
```@example mixing
d = mixing(1.0, 0.09, 1.25e-3, 2.2, 2e3, 0.04, 1100.0, 0.001, 0.0, 0.0)
c + d
```
Given a vector or other collection of samples, they can be combined by summing them.
```@example mixing
v = [a, b, c, d]
sum(v)
```
These operations are also fast because they are just doing arithmetic and forming named tuples underneath.
```@example mixing
@benchmark sum($v)
```

## A Simple Simulation

The mixing function itself can be used for simulations, instead of an entire hypothetical deployment as in the [tutorial](tutorial.md) and [efficient tutorial](efficiency.md). For example, we can use the [Turing.jl](https://turing.ml/) modeling tools to do a simple Monte Carlo simulation, looking at how parameters in the mixing model relate to each other.

Here I define a model where each parameter is assigned a probability distribution to sample from. To make it as simple as possible, I just use uniform distributions for each parameter, drawing randomly from a prescribed range.

```@example mixing
using Turing

import CairoMakie as Mke

@model function model()

    γ ~ Uniform()
    d ~ Uniform(0.05, 0.15)
    a = 1.25e-3
    Q ~ Uniform(1.0, 4.0)
    ρf = 2e3
    cf ~ Uniform(4e4, 6e4) # defined in ppm for convenience
    ρs ~ Uniform(750, 1250)
    cs ~ Uniform(1e3, 4e3) # defined in ppm for convenience
    𝓁 ~ Uniform()
    ℒ ~ Uniform()

    core = mixing(γ, d, a, Q, ρf, cf / 1e6, ρs, cs / 1e6, 𝓁, ℒ)
    # capture the core's mixed concentration
    c ~ Dirac(1e6 * core[:analyte])

end
```
This model can be sampled like any model. In this case, we use the prior sampler to do the Monte Carlo for us. Below, we create a chain of 64,000 samples and then print out the summary.
```@example mixing
chain = sample(model(), Prior(), 64_000)
chain # hide
```
We can do anything we want with these samples. Here are a few plots showing the relationships between each parameter and the core concentration after mixing.
```@example mixing
fig = Mke.Figure(size=(700,1100))
for (i,x) in enumerate([:γ, :d, :Q, :cf, :ρs, :cs, :𝓁, :ℒ])
    ax = Mke.Axis(
        fig[(i - 1) ÷ 2 + 1,
        (i - 1) % 2 + 1],
        xlabel=string(x),
        ylabel="core concentration (ppm)"
    )
    Mke.hexbin!(ax, chain[x] |> vec, chain[:c] |> vec, bins=30)
end
fig
```
There's one very clear feature of the plots above: the baseline soil concentration (`cs`) absolutely dominates the response in the core concentration. There is an almost 1:1 correspondence between baseline soil concentration and mixed soil concentration. To see the influence of other parameters, we could hold baseline soil constant and vary everything else around it.
