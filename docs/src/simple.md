# Simple Simulation & Analysis

This page is an example of a simple simulation.

By "simple" I mean that we ignore any potential spatial covariance in any of the mixing and deployment parameters. The specific location of each soil core/sample, therefore, doesn't matter. We treat them as fully independent realizations of the same mixing process. To analyze the simulation, we look at the difference in average elemental concentration between two different sampling rounds. In doing so, we get a sense of how reliably we can resolve weathering and tracer signals with different numbers of samples.

## Setup

```@example simple
using Monty, Meshes, Random, Distributions, DataFrames

import CairoMakie as Mke

const rng = Xoshiro(6)
nothing # hide
```

Because we forgo spatial relationships among sample locations, we don't care about their locations at all. They can just be drawn from any arbitrary geometry. To make it simple, the deployment geometry can just be a box. It doesn't matter.
```@example simple
geometry = Box((0,0), (1,1))
geometry # hide
```
We'll take two rounds of samples. The first round is exactly after spreading and the next round is 8 months later.
```@example simple
times = [0.0, 8.0/12.0]
nothing # hide
```
Assume the sampling depth varies between 10 and 20 cm. Later, we'll also assume that the feedstock isn't mixed deeper than 10 cm.
```@example simple
depth = TriangularDist(0.1, 0.2)
nothing # hide
```
Now we define distributions that mixing parameters are drawn from for each core. We'll use two generic analytes, one called "cation" and the other called "tracer." The background soil concentrations for the cation and tracer vary about 20 % ($\sigma/\mu$), as defined by the Gamma distributions below.
```@example simple
# application rate [kg/m^2]
Q = Gamma(20.0, 0.1)
# feedstock bulk density [kg/m^3]
ρf = 2000.0
# feedstock elemental concentrations [kg/kg]
cf = (
    cation=Normal(0.05, 0.05/10),
    tracer=Normal(0.001, 0.001/10)
)
# soil bulk density [kg/m^3]
ρs = Normal(900, 100)
# soil elemental concentrations [kg/kg]
cs = (
    cation=Gamma(20, 0.0005),
    tracer=Gamma(20, 0.0000125)
)
nothing # hide
```
For the leaching models, I assume simple exponential leaching for the cation and no leaching for the tracer. For an exponential decay parameter of 1.0, we expect just about half of the cation to be lost over 8 months.
```@example simple
leach = (
    cation=ExponentialLeaching(λ=1.0),
    tracer=NoLeaching()
)
nothing # hide
```

# Simulation

We'll create simulation stacks with different sample sizes, to see how the number of samples influences an estimate of the difference in mean concentrations across the two sampling rounds. First, I'll run through a single simulation.

Sample locations are random in each sampling round, but again, the sample locations do not matter when there is no spatial correlation or other spatial influence (think Monty version 1). I'll use 6 cores per composite sample and 100 total samples. I also use 4 % and 7 % measurement error for the cation and tracer, respectively.
```@example simple
plan = randomsampleplan(rng, geometry, 100, times)
samp = CoreSet(plan, 6)
sim = Simulation((:cation, :tracer), samp)
nothing # hide
```
Then we just fill in the fields.
```@example simple
unmixed!(rng, sim, depth=depth)
rand!(rng, Q, sim.Q)
sim.ρf .= ρf
rand!(rng, ρs, sim.ρs)
for analyte ∈ (:cation, :tracer)
    rand!(rng, cf[analyte], sim.cf[analyte])
    rand!(rng, cs[analyte], sim.cs[analyte])
    leaching!(sim, analyte, leach[analyte], plan.time)
end
massloss!(sim, plan, mean)
# compositing and measurement error
analyze!(rng, sim, (cation=0.04, tracer=0.07), 0.005)
df = GeoTable(sim, plan) |> DataFrame
describe(df)
```

## Repeated Simulation

Now we run repeated simulations with different numbers of samples by wrapping the simulation script in the [`simulationstack`](@ref) function. In each case, we still use 6 cores per composite and do 10,000 simulations per case. Inside this loop, I also convert the results for each stack into a DataFrame for more familiar processing. Each of these simulation stacks are finished (on my computer) in under one second.

```@example simple
sample_sizes = [20, 40, 60, 80, 100]
stacks = Vector{Any}(undef, length(sample_sizes))
for (i,n) ∈ sample_sizes |> enumerate
    plan = randomsampleplan(rng, geometry, n, times)
    samp = CoreSet(plan, 6)
    sim = Simulation((:cation, :tracer), samp)
    sims = simulationstack(10_000, sim, samp, plan, show_progress=false) do
        unmixed!(rng, sim, depth=depth)
        rand!(rng, Q, sim.Q)
        sim.ρf .= ρf
        rand!(rng, ρs, sim.ρs)
        for analyte ∈ (:cation, :tracer)
            rand!(rng, cf[analyte], sim.cf[analyte])
            rand!(rng, cs[analyte], sim.cs[analyte])
            leaching!(sim, analyte, leach[analyte], plan.time)
        end
        massloss!(sim, plan, mean)
        analyze!(rng, sim, (cation=0.05, tracer=0.07), 0.005)
    end
    stacks[i] = sims
end
nothing # hide
```


## Analysis

To process and analyze these results, I convert each simulation stack to a DataFrame and get rid of unnecessary columns. Then I group each frame and average the analyte concentrations in each round before taking their difference. The result is a group of new data frames with the *difference in mean* values for each analyte between sampling rounds 1 and 2, for all realizations. I know this is a dense little code block, but the focus here is not dataframe functionality, which is documented [here](https://dataframes.juliadata.org/stable/man/split_apply_combine/) in [DataFrames.jl](https://dataframes.juliadata.org/stable/).

```@example simple
difs = map(stacks) do sims
    combine(
        groupby(
            select(sims |> DataFrame, :realization, :analyte, :round, :data),
            [:realization, :analyte],
        ),
        [:round, :data] => (r, d) -> mean(d[r.==2]) - mean(d[r.==1]),
    )
end
nothing # hide
```

We can plot these differences in means for each sample size (I hide the code and just show the plot). These densities show us the relative likelihood of differences difference in means for the mobile cation, the tracer, and the sample mass, given all of our input parameters.
* We expect the cation concentration to go *down* over time, between sampling rounds 1 and 2, because it's weathering and leaching. So, these densities should be negative, on average. They are, and as the sample size increases, the width of the density shrinks.
* The tracer concentration change should be about zero because we assumed it doesn't leach out of the soil at all. The tracer densities look fine, all with their mean values at approximately zero.
* We don't care as much about sample masses, but those should also be approximately the same for each sampling round, so the difference should be zero. The densities look fine.
These kinds of comparisons can be used to estimate the number of samples required for a certain level of confidence in average concentration changes between sampling rounds. They do not, quite yet, tell you how many samples are required for a desired level of confidence in initial CDR. For that, you have to apply your method of estimating/inferring CDR to these datasets and test the sensitivity.
```@example simple
fig = Mke.Figure(size=(1200, 1200)) # hide
for i ∈ eachindex(sample_sizes) # hide
    for (j, (analyte, g)) ∈ pairs(groupby(difs[i], :analyte)) |> enumerate # hide
        s = analyte[:analyte] # hide
        x = g[!, :round_data_function] # hide
        title = # hide
            s == :mass ? "$s concentration change [kg]\n$(sample_sizes[i]) samples" : # hide
            "$s concentration change [ppm]\n$(sample_sizes[i]) samples" # hide
        if j == 1 # hide
            limits = ((-2000, 2000), nothing) # hide
        elseif j == 2 # hide
            limits = ((-50, 50), nothing) # hide
        else # hide
            limits = ((-0.1, 0.1), nothing) # hide
        end # hide
        ax = Mke.Axis(fig[i, j], title=title, limits=limits) # hide
        Mke.density!(ax, s == :mass ? x : 1e6 * x) # hide
        Mke.vlines!(ax, 0, color=:black) # hide
    end # hide
end # hide
fig # hide
```
