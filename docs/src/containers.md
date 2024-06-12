# Containers

`Monty` simulations operate on three main data types defined by the package, which function as containers for different information. These container types are useful because they abstract away some details when running simulations and because they only allow certain operations to be done on them, which helps enforce correctness. They are described in the sections below in detail. Doing very specific simulations or modifying the package will require some understanding of these container types and their internals.

## `SamplePlan`

The [`SamplePlan`](@ref) defines the ideal outcome for a sampling strategy. Its docstring is presented below.

```@docs
SamplePlan
```

As the docstring indicates, a sample plan contains fields for the location, sampling round, time, control status, and coordinates (points) of each point in the sampling plan. These fields are all flat vectors (1-dimensional) and a `SamplePlan` functions is much like a table with these vectors in its columns. Because of this, it can easily be converted to a `GeoTable`. Also, note that the fields of a sample plan are `ReadOnly`. They can't be modified once constructed.

To explain these vectors more fully:
* The `location` field is an integer index that indicates the location status of a sample. This is usually only relevant for paired sample plans when locations are supposed to be revisited. For example, a 2 x 2 grid of points has 4 locations, and each of these locations will be repeated in each sampling round. If there were two sampling rounds, the `location` vector might look like: `[1, 2, 3, 4, 1, 2, 3, 4]`. This means that each location is visited twice, one time in each sampling round.
* The `round` vector is also an integer index. It indicates the ordering/timing of the sample. For example, if there were two rounds of sampling on a 2 x 2 grid of points, the `round` vector would look like `[1, 1, 1, 1, 2, 2, 2, 2]`.
* The `time` vector indicates the timing of each sample. Samples in the same `round` will have identical times and there is a one-to-one relationship between `round` and `time`.
* The `control` vector contains booleans (`true`/`false`) indicating whether samples are controls. Controls do not get feedstock applied to them after `time = 0`.
* The `points` vector contains the two-dimensional coordinates of each point in the plan. These are [`Point`](https://juliageometry.github.io/MeshesDocs/stable/geometries/primitives.html#Meshes.Point) types from the [`Meshes.jl`](https://juliageometry.github.io/MeshesDocs/stable/index.html).

A `SamplePlan` can be created by defining all of these fields directly, but this is not generally an efficient way to do it. For example, here is a very simple plan for two cells of a grid, one of which is a control cell
```@example containers
using Monty
using Meshes

# locations are repeated across two rounds of sampling
location = UInt16[1, 2, 1, 2]
# there are two sampling rounds for two locations
round = UInt16[1, 1, 2, 2]
# the times are arbitrary, but there are only two of them for the two sampling rounds
time = [0.0, 0.0, 1.0, 1.0]
# one cell is control, the other is treatment, across both sampling rounds
control = [true, false, true, false]
# the two cell's coordinates are also repeated
points = [Point(1.0, 1.0), Point(1.0, 2.0), Point(1.0, 1.0), Point(1.0, 2.0)]

SamplePlan(location, round, time, control, points)
```

It's usually easier to use the sample planning functions, [`pairedsampleplan`](@ref) and [`randomsampleplan`](@ref), to generate a valid sample plan. For example, using a 2 x 2 grid
```@example containers
grid = CartesianGrid(2, 2)
plan = pairedsampleplan(grid .|> centroid, [0.0, 1.0])
plan # hide
```

## `CoreSet`

The [`CoreSet`](@ref) is the second container needed for simulations and it's very simple. It's just a wrapper around a two-dimensional array (a matrix) of points. Its definition is just this:
```
struct CoreSet{𝒯}
    points::Matrix{Point{2,𝒯}}
end
```
The `points` field contains the coordinates of *each core* when a sample plan is executed. Each column of the matrix contains the core locations that will be composited into an individual sample. There should be exactly as many columns of the `points` matrix as there are samples in the `SamplePlan` for a simulation.

Generally, it should be easiest to allocate a `CoreSet` by calling the constructor on an existing `SamplePlan`. For example, using the one from above and specifying 10 cores per sample:
```@example containers
coreset = CoreSet(plan, 10)
coreset # hide
```
Initially, the core set is just empty. This is what the internal matrix looks like:
```@example containers
coreset.points
```

## `Simulation`

The [`Simulation`](@ref) container is the big one. Here is the docstring:
```@docs
Simulation
```
The `nsamp` and `ncore` fields are just integers indicating the total number of samples and cores in the plan. The `a` field is a scalar indicating the cross-sectional area of the soil cores. Most of the other fields are matrices (two-dimensional arrays) of numbers *or* named tuples of matrices. If the fields apply to multiple analytes, they contain named tuples of matrices where the keys/names indicated the analyte names. The `cores` field is a matrix of soil [`Sample`](@ref) types and the `composites` field is a flat vector for the results of compositing/combining the `cores`. The `measurements` field is also a vectors, which stores the results of applying random measurement error to the `composites`.

Allocating a simulation container is usually easiest with a pre-existing `CoreSet`. For example, using the one from above:
```@example containers
sim = Simulation((:Ca, :Mg, :Fe, :Cr), coreset)
sim # hide
```
The printed summary for a `Simulation` is helpful, but intentionally obscures the structure of internal fields. Any of the fields can be accessed directly, however. Here is what the `cs` field, which is a named tuple of matrices for soil concentrations, looks like for our `:Ca` analyte (calcium). There are 10 rows for the 10 cores in each composite sample and 8 columns for 8 total composite samples:
```@example containers
sim.cs[:Ca]
```
This is what the `𝓁` field (loss fractions) looks like for the `:Fe` analyte, also just a matrix full of `NaN`s:
```@example containers
sim.𝓁[:Fe]
```
and so on...

All of these fields are initially empty or `NaN`. Performing a simulation requries filling in all of these fields, `γ` through `ℒ`. They represent the mixing and weathering parameters for the [mixing model](mixing.md). There are helper functions to fill some of these fields in, like [`spreading!`](@ref), [`soilconcentration!`](@ref), [`leaching!`](@ref), mixing functions like [`triangularmixing!`](@ref), and others. Once all the fields are filled, we call the [`analyze!`](@ref) function to automatically fill in `cores`, `composites`, and `measurements`

What matters is really *how* the mixing fields are filled in and this is up to you. They can be filled with a single value everywhere, by drawing from a distribution, by drawing from a [`GaussianSimulator`](@ref) or [`GaussianCosimulator`](@ref), or really by any method that aligns with some kind of simulation goal.

The [simulation tutorial](tutorial.md) page goes through one example in detail.
