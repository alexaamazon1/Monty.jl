# Monty

### Overview

Monty is a package for simulating geochemical data collected for the quantification of carbon dioxide removal (CDR) in enhanced rock weathering (ERW) field trials and commercial deployments. The package simulates whole-rock chemical composition measurements used in mass-balance methods. It simulates *both* the data collection process *and* physical processes/properties *jointly*, to produce realistic synthetic datasets that account for interactions between operational choices and physical/chemical conditions. It's a sandbox for mass-balance MRV in ERW.

### Mixing Model

The core element of the package is a simple deterministic *mixing model* of how the chemical concentration of an *individual soil core* is determined by
1. initial/background soil chemical composition
2. soil bulk density
3. feedstock chemical composition
4. feedstock application rate
5. feedstock mixing profile with depth
6. feedstock bulk density
7. dissolution and leaching of mobile elements over time
8. total mass loss from feedstock
9. soil core depth

This mixing model yields a [composite type](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) that contains the *exact* (more on this later) concentrations of each analyte in the collected soil mass and the *exact* mass of the collected soil. Every evaluation represents a single, perfect soil core.

Look at the [mixing tutorial page](mixing.md) for examples of running the mixing model.

### Compositing

*Sample compositing*—the common practice of physically combining some number of soil cores into one sample—is performed using the addition function (`+`), making it trivially easy to aggregate any number of cores/samples using `+` or any function that performs addition like `sum`.

### Sample Planning

Aside this central mixing model is a group of modules for simulating the planning and collection of soil samples in time and space. These modules help define field/deployment geometry and the time/location of soil cores taken from hypothetical deployments. They handle processes like
* selecting sample points from treatment & control geometries
  * cartesian grids or arbitrary polygons
  * repeating locations for each sampling round (paired sampling) or not (independent sampling)
* assigning times and control labels to samples
* simulating the error/jitter of real sampling, where the prescribed locations are only approximately visited
* defining the pattern/stencil used for core compositing (physical combination) at each sample location

### Simulation

Finally, on top of both these components, are functions and types for simulating data produced by the physical/chemical parameters above and sampling choices *with prescribed variability*. By treating all the inputs as stochastic parameters (drawn from a known probability distribution) simulations capture the effects of variability in all these physical and operational parameters at the same time.

### Measurement

At the end of a simulation, composited soil samples are "measured" with analytical error, usually defined in terms of a relative standard deviation, *closely representing the real data collection process from beginning to end*. As mentioned above, sample concentrations may vary according to prescribed parameters and distributions, but until the measurement step they are represented *exactly*. Applying random measurement noise obscures these concentrations. The language can be fuzzy, but the pedantic division would be between physical and chemical *variability* in the sample population and *noise* in our measurements.

### Implementation

The primary concern for the code itself is correctness—encouraging and enforcing it, minimizing opportunities for mistakes. Physical quantities defined for elements/analytes are generally done so with [named tuples](https://docs.julialang.org/en/v1/manual/types/#Named-Tuple-Types). This allows various functions to require the chosen group of analytes *and their ordering* to be consistent. Named tuples are also fast in this context, avoiding memory allocation. Additionally, almost all functions are strict about the use of consistent data types and liberal with checks and assertions.

# Why simulate?

A few reasons:
1. **Understanding and organization**. Simulating a simple *forward model* of a physical/chemical process that we're trying to measure in the field illuminates what information is required to measure the process accurately and precisely. By putting together a model of the *data generating-process* we're forced to think about which parameters are relevant, how they interact, and what we expect to see in real data. Writing organizes your thoughts.
2. **Uncertainty analysis**. Simulations directly address questions about how variability in the modeled system translates to uncertainty in statistical results. We can simulate hypothetical deployments with different levels of variability and see how it impacts our confidence in statistical results based on the simulated data.
3. **Statistical validation**. The forward model can be used to generate synthetic datasets for which *we know the right answer unambiguously*. These synthetic datasets can be used to *check* the statistical methods we use for real data by running those methods on the synthetic data and seeing if they fail. For example, *before doing the trial*, we could follow roughly these steps to validate a statistical method for our sampling plan and soil/feedstock properties:
    - define Monty parameters that represent a hypothetical deployment with $X$ tons of CDR, on average
    - generate 10,000 realizations of this deployment with prescribed levels of uncertainty and variability
    - apply the statistical method of choice to every synthetic dataset for this hypothetical deployment
    - check that CDR estimates are, on average, very close to $X$ tons and, very importantly, check that the uncertainty quantification associated with these calculations *is consistent with expectations*
4. **Power analysis**. Synthetic data can also be used to appropriately *power* the sampling plan for field trials or deployments. By making conservative assumptions about variability and noise (or collecting appropriate reconnaissance data) and simulating with those assumptions, sample sizes, and sampling plans can be selected based on the uncertainty associated with synthetic data. Put more plainly, we can do the following
    - make assumptions about noise and variability
    - simulate with different sample sizes
    - analyze the simulated data as if it were the real data collected from the experiment/trial
    - see what sample size is appropriate for a target level of uncertainty

# Workflow

Simulating involves a few basic steps, each of which can be approached differently but all of which should be interoperable/composable. For a complete example, see [the example](tutorial.md). For a demo of efficient repeated examples, see the [efficient simulations](efficiency.md) page.

1. A hypothetical deployment is defined in terms of a sample plan. This represents the exact location and time of every sample that will be taken, including whether they are control samples that will not have feedstock spread on their locations. This can be done with almost any field geometry and a number of common sampling layouts, in addition to random sample locations. There are several functions for creating sample plans with *paired* or *independent* locations in each round of sampling. Feedstock application is assumed to occur exactly at time zero.
2. The sample plan is *realized*. This means sample locations are visited by an imaginary sampler, who isn't perfectly accurate, and each sample is broken out into some number of core locations that will eventually be composited into a single sample.
3. Physical/chemical parameters for *each core* in the deployment are defined, probably by drawing them from prescribed distributions. Properties like baseline soil concentrations and application rates can easily be defined with spatial autocorrelation.
4. Simulated properties for each core are passed to the mixing model, generating sample concentrations. Then the individual cores for each location and time are composited and measured, producing the information that would ultimately be observed for the hypothetical deployment.

# Reference

See the [reference page](reference.md) for documentation of types and functions.
