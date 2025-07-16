module Sampling

using Random: AbstractRNG, randn, rand
using Statistics: mean
using Unitful: @u_str
using UnPack: @unpack
using ReadOnlyArrays: ReadOnlyArray
using Meshes:
    Geometry,
    GeometrySet,
    Mesh,
    Point,
    PointSet,
    Vec,
    Box,
    CartesianGrid,
    coordinates,
    centroid
using GeoTables: georef

import GeoTables: GeoTable
import Meshes: boundingbox
import Base: show

using ..Stencils
using ..Jitters

import ..ncore, ..nsample

using DocStringExtensions

#--------------------------------------
export randin, randin!

"""$(TYPEDSIGNATURES)
Generates a random point inside a two-dimensional geometry by sampling randomly in the bounding box until the point is contained by the geometry. No more than `maxiter` samples are drawn."""

function randin(
    rng::AbstractRNG,
    geom::Geometry{2,𝒯};
    maxiter::Integer=1_000_000,
)::Point{2,𝒯} where {𝒯}
    box = boundingbox(geom)
    x₁, y₁ = box |> minimum |> coordinates
    x₂, y₂ = box |> maximum |> coordinates
    Δx, Δy = (x₂ - x₁), (y₂ - y₁)
    point = Point(x₁ - one(𝒯), y₁ - one(𝒯)) # initial point outside the geometry
    iter::Int64 = 0
    while !(point ∈ geom)
        x = Δx * rand(rng, 𝒯) + x₁
        y = Δy * rand(rng, 𝒯) + y₁
        point = Point(x, y)
        iter += 1
        if iter > maxiter
            error(
                "Point sampling from inside the provided geometry ($(typeof(geom))) reached $maxiter iterations. Make sure you used a geometry that represents 2D area and not just a line (e.g. a Ball{2} will work but not a Sphere{2}).",
            )
        end
    end
    return point
end

"""$(TYPEDSIGNATURES)
Fills a vector with random points inside a two-dimensional geometry"""
function randin!(
    rng::AbstractRNG,
    points::AbstractVector{Point{2,𝒯}},
    geom::Geometry{2,𝒯},
)::Nothing where {𝒯}
    for i ∈ eachindex(points)
        points[i] = randin(rng, geom)
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Generates a vector of `N` random points inside a two-dimensional geometry"""
function randin(
    rng::AbstractRNG,
    geom::Geometry{2,𝒯},
    N::Int,
)::Vector{Point{2,𝒯}} where {𝒯}
    points = Vector{Point{2,𝒯}}(undef, N)
    randin!(rng, points, geom)
    return points
end

#--------------------------------------
export gridcelloverlay, gridpointoverlay

function gridcelloverlay(geom, dims::NTuple{2,Integer})
    box = geom |> boundingbox
    grid = CartesianGrid(box |> minimum, box |> maximum, dims=dims)
    filter(!isnothing, (grid |> collect) .∩ geom)
end

"""$(TYPEDSIGNATURES)
Creates a [`CartesianGrid`](https://juliaearth.github.io/GeoStatsDocs/stable/domains.html#CartesianGrid) of size `dims`, intersects the grid with `geom`, and returns the centroids of intersected grid cells. Any cells with no intersection are removed."""
function gridpointoverlay(geom, dims::NTuple{2,Integer})
    gridcelloverlay(geom, dims) .|> centroid
end

"""$(TYPEDSIGNATURES)
Creates a [`CartesianGrid`](https://juliaearth.github.io/GeoStatsDocs/stable/domains.html#CartesianGrid) of size `dims`, intersects the grid with `geom`, and returns points *sampled randomly* from the intersected grid cells. Any cells with no intersection are removed."""
function gridpointoverlay(rng::AbstractRNG, geom, dims::NTuple{2,Integer})
    [randin(rng, cell) for cell ∈ gridcelloverlay(geom, dims)]
end

#------------------------------------------------------------------------------
export SamplePlan, GeoTable, boundingbox

_sample_plan_input(x, sortidx) = view(x, sortidx) |> copy |> ReadOnlyArray

"""A `SamplePlan` represents the intended locations and times for all samples in a simulated deployment, in addition to whether target samples are control samples (no feedstock spreading). Generally, it's easiest to construct a sample plan using either the [`pairedsampleplan`](@ref) or [`randomsampleplan`](@ref) functions. The fields of a `SamplePlan` are [`ReadOnlyArray`](https://github.com/bkamins/ReadOnlyArrays.jl) types that cannot be modified.

The fields are
$(TYPEDFIELDS)
"""
struct SamplePlan{𝒯}
    location::ReadOnlyArray{UInt16,1,Vector{UInt16}}
    round::ReadOnlyArray{UInt16,1,Vector{UInt16}}
    time::ReadOnlyArray{𝒯,1,Vector{𝒯}}
    control::ReadOnlyArray{Bool,1,Vector{Bool}}
    points::ReadOnlyArray{Point{2,𝒯},1,Vector{Point{2,𝒯}}}
    function SamplePlan(
        location::Vector{UInt16},
        round::Vector{UInt16},
        time::Vector{𝒯},
        control::Vector{Bool},
        points::Vector{Point{2,𝒯}},
    ) where {𝒯}
        @assert length(location) ==
                length(round) ==
                length(time) ==
                length(control) ==
                length(points)
        @assert minimum(round) == 1
        @assert all(i -> i ∈ round, 1:maximum(round))
        @assert minimum(location) == 1
        @assert all(i -> i ∈ location, 1:maximum(location))
        @assert issorted(time)
        # the location indices will be sorted
        idx = sortperm(location)
        new{𝒯}(
            _sample_plan_input(location, idx),
            _sample_plan_input(round, idx),
            _sample_plan_input(time, idx),
            _sample_plan_input(control, idx),
            _sample_plan_input(points, idx),
        )
    end
end

"""$(TYPEDSIGNATURES)
Returns the number of sample locations in a [`SamplePlan`](@ref). This is the number of planned samples *after compositing* the cores for each planned location. The number of cores *for each sample* is given by `ncore(plan)`. The total number of cores is `ncore(plan) * nsample(plan)`."""
nsample(plan::SamplePlan) = length(plan.points)

"""$(TYPEDSIGNATURES)
Returns the bounding box of the samples in a `plan`"""
boundingbox(plan::SamplePlan) = boundingbox(plan.points |> PointSet)

function show(io::IO, plan::SamplePlan{𝒯}) where {𝒯}
    @unpack location, round, time, points = plan
    println(io, "SamplePlan{$𝒯)")
    println(io, "├─ total samples = $(nsample(plan))")
    println(io, "├─ sampling rounds = $(maximum(round))")
    println(io, "├─ sampling times = $(sort(unique(time)))")
    println(io, "└─ bounding box")
    println(io, "  ├─ min = $(points |> PointSet |> boundingbox |> minimum)")
    println(io, "  └─ max = $(points |> PointSet |> boundingbox |> maximum)")
    nothing
end

"""$(TYPEDSIGNATURES)
Creates a tabular representation of the sample plan using the [`GeoTable`](https://juliaearth.github.io/GeoStatsDocs/stable/data.html) format."""
function GeoTable(plan::SamplePlan)
    @unpack location, round, time, control, points = plan
    georef((location=location, round=round, time=time, control=control), points)
end

#--------------------------------------
export pairedsampleplan

function fillpairedplan(points, times)
    location = repeat(collect(UInt16, 1:length(points)), length(times))
    rounds = repeat(collect(UInt16, 1:length(times)), inner=length(points))
    points, times =
        repeat(points, length(times)), repeat(times, inner=length(points))
    return location, rounds, points, times
end

"""$(TYPEDSIGNATURES)
Creates a [`SamplePlan`](@ref) where sample locations are repeated for each round of sampling (location-paired). No control points are created. The timing of each sampling round is given by the `times` vector."""
function pairedsampleplan(
    points::Vector{Point{2,𝒯}},
    times::AbstractVector{𝒯},
)::SamplePlan{𝒯} where {𝒯}
    location, rounds, points, times = fillpairedplan(points, times)
    controls = fill(false, length(points))
    SamplePlan(location, rounds, times, controls, points)
end

"""$(TYPEDSIGNATURES)
Creates a location-paired smaple plan using a vector of existing `points`. Points inside the `control` geometry will be assigned to the control group. The timing of each sampling round is given by the `times` vector."""
function pairedsampleplan(
    points::Vector{Point{2,𝒯}},
    control,
    times::AbstractVector{𝒯},
)::SamplePlan{𝒯} where {𝒯}
    location, rounds, points, times = fillpairedplan(points, times)
    controls = map(point -> point ∈ control, points)
    SamplePlan(location, rounds, times, controls, points)
end

"""$(TYPEDSIGNATURES)
Creates a location-paired sample plan using existing vectors of `treatment` and `control` points. The timing of each sampling round is given by the `times` vector."""
function pairedsampleplan(
    treatments::Vector{Point{2,𝒯}},
    controls::Vector{Point{2,𝒯}},
    times::AbstractVector{𝒯},
)::SamplePlan{𝒯} where {𝒯}
    L = length(treatments) + length(controls)
    ntime = length(times)
    SamplePlan(
        repeat(collect(UInt16, 1:L), ntime),
        repeat(collect(UInt16, 1:ntime), inner=L),
        repeat(times, inner=L),
        repeat(
            [zeros(Bool, length(treatments)); ones(Bool, length(controls))],
            ntime,
        ),
        repeat([treatments; controls], ntime),
    )
end

"""$(TYPEDSIGNATURES)
Creates a location-paired sample plan using `N` points selected random from inside a `field` geometry. The locations are repeated for each sampling round specified in the `times` vector."""
function pairedsampleplan(
    rng::AbstractRNG,
    field::Geometry{2,𝒯},
    N::Integer,
    times::AbstractVector{𝒯},
) where {𝒯}
    pairedsampleplan(randin(rng, field, N), times)
end

"""$(TYPEDSIGNATURES)
Creates a [`SamplePlan`](@ref) where sample locations are repeated for each round of sampling (location-paired). `ntreatment` treatment locations are drawn randomly from inside the `treatment` geometry and `ncontrol` control locations are drawn randomly from inside the `control` geometry. The timing of each sampling round is given by the `times` vector."""
function pairedsampleplan(
    rng::AbstractRNG,
    treatment::Geometry{2,𝒯},
    ntreatment::Integer,
    control::Geometry{2,𝒯},
    ncontrol::Integer,
    times::AbstractVector{𝒯},
) where {𝒯}
    pairedsampleplan(
        vcat(
            randin(rng, treatment, ntreatment),
            randin(rng, control, ncontrol),
        ),
        control,
        times,
    )
end

#--------------------------------------
export randomsampleplan

"""$(TYPEDSIGNATURES)
Creates a random [`SamplePlan`](@ref), where sample locations are not repeated across different rounds of sampling. For each round of sampling, treament points are drawn randomly from the `treatment` geometry. The number of treatment points in each round is determined by the vector of integers `ntreatment`. The same true for control points, using the `control` geometry and `ncontrol` sample sizes for each round. The timing of each sampling round is given by the `times` vector."""
function randomsampleplan(
    rng::AbstractRNG,
    treatment::Geometry{2,𝒯},
    ntreatment::AbstractVector{<:Integer},
    control::Geometry{2,𝒯},
    ncontrol::AbstractVector{<:Integer},
    times::AbstractVector{𝒯},
) where {𝒯}
    @assert length(ntreatment) == length(ncontrol) == length(times)
    # number of sampling rounds
    nround = length(ntreatment)
    # total number of target locations
    L = sum(ntreatment) + sum(ncontrol)
    # empty vectors
    rounds = Vector{UInt16}(undef, L)
    points = Vector{Point{2,𝒯}}(undef, L)
    controls = Vector{Bool}(undef, L)
    a = 1
    for i ∈ 1:nround
        # fill treatment vectors for round i
        b = a + ntreatment[i] - 1
        rounds[a:b] .= i
        randin!(rng, view(points, a:b), treatment)
        controls[a:b] .= false
        # fill control vectors for round i
        a += ntreatment[i]
        b = a + ncontrol[i] - 1
        rounds[a:b] .= i
        randin!(rng, view(points, a:b), control)
        controls[a:b] .= true
        # increment for next round
        a += ncontrol[i]
    end
    SamplePlan(collect(UInt16, 1:L), rounds, times[rounds], controls, points)
end

"""$(TYPEDSIGNATURES)
Creates a random [`SamplePlan`](@ref), where sample locations are not repeated across different rounds of sampling. For each round of sampling, `ntreatment` treament points are drawn randomly from the `treatment` geometry and `ncontrol` points are drawn randomly from inside the `control` geometry. The timing of each sampling round is given by the `times` vector."""
function randomsampleplan(
    rng::AbstractRNG,
    treatment::Geometry{2,𝒯},
    ntreatment::Integer,
    control::Geometry{2,𝒯},
    ncontrol::Integer,
    times::AbstractVector{𝒯},
) where {𝒯}
    randomsampleplan(
        rng,
        treatment,
        fill(ntreatment, length(times)),
        control,
        fill(ncontrol, length(times)),
        times,
    )
end

"""$(TYPEDSIGNATURES)
Creates a random sample plan with points selected randomly from inside the `field` geometry for each round. The number of points in each round is given by a vector of integers `N` and the timing of each round by the `times` vector."""
function randomsampleplan(
    rng::AbstractRNG,
    field::Geometry{2,𝒯},
    N::AbstractVector{<:Integer},
    times::AbstractVector{𝒯},
) where {𝒯}
    @assert length(N) == length(times)
    L = sum(N)
    rounds = Vector{UInt16}(undef, L)
    points = Vector{Point{2,𝒯}}(undef, L)
    controls = Vector{Bool}(undef, L)
    a = 1
    for i ∈ 1:length(times)
        # fill vectors for round i
        b = a + N[i] - 1
        rounds[a:b] .= i
        randin!(rng, view(points, a:b), field)
        controls[a:b] .= false
        # increment for next round
        a += N[i]
    end
    SamplePlan(
        collect(UInt16, 1:L),
        rounds,
        times[rounds],
        fill(false, L),
        points,
    )
end

"""$(TYPEDSIGNATURES)
Creates a random [`SamplePlan`](@ref), where sample locations are not repeated across different rounds of sampling. For each round of sampling, `N` treament points are drawn randomly from the `field` geometry. The timing of each sampling round is given by the `times` vector."""
function randomsampleplan(
    rng::AbstractRNG,
    field::Geometry{2,𝒯},
    N::Integer,
    times::AbstractVector{𝒯},
) where {𝒯}
    L = length(times) * N
    rounds = repeat(collect(UInt16, 1:length(times)), inner=N)
    times = repeat(times, inner=N)
    controls = fill(false, L)
    points = randin(rng, field, L)
    SamplePlan(collect(UInt16, 1:L), rounds, times, controls, points)
end

#--------------------------------------
export CoreSet, executeplan, executeplan!

"""A `CoreSet` represents the *realization* of a [`SamplePlan`](@ref) and is closely related. The sample plan represents the ideal sampling locations. A core set expands on the sample plan by
1. representing planned and unplanned jitter in the target sample locations (for example, simple GPS inaccuracy)
2. defining some number of *cores* for each sample, which are composited during the simulation process.

Dividing the sample plan and core set allows for more efficient simulation of many realizations from the same hypothetical sample plan. A core set can be defined without any location jitter and without compositing (a single core).

A `CoreSet` can be created and filled by the [`executeplan`](@ref) and [`executeplan!`](@ref) functions, respectively."""
struct CoreSet{𝒯}
    points::Matrix{Point{2,𝒯}}
end

function CoreSet(T::Type, samples::Integer, cores::Integer)
    CoreSet(fill(Point(zero(T), zero(T)), cores, samples))
end

"""$(TYPEDSIGNATURES)
Allocates a `CoreSet` with the correct size by referencing a sample `plan` with `n` cores per sample."""
function CoreSet(plan::SamplePlan{𝒯}, cores::Integer) where {𝒯}
    CoreSet(𝒯, plan |> nsample, cores)
end

"""$(TYPEDSIGNATURES)
Allocates a `CoreSet` with the correct size by referencing a sample `plan` for the number of samples and a `stencil` for the number of cores at each sample."""
function CoreSet(plan::SamplePlan{𝒯}, stencil::AbstractStencil{𝒯}) where {𝒯}
    CoreSet(𝒯, plan |> nsample, stencil |> ncore)
end

"""$(TYPEDSIGNATURES)
Returns the total number of sample locations in a `CoreSet`"""
ncore(samp::CoreSet) = size(samp.points, 1)

"""$(TYPEDSIGNATURES)
Returns the number of cores per sample"""
nsample(samp::CoreSet) = size(samp.points, 2)

boundingbox(samp::CoreSet) = boundingbox(samp.points |> vec |> PointSet)

function show(io::IO, samp::CoreSet{𝒯}) where {𝒯}
    @unpack points = samp
    println(io, "CoreSet{$𝒯}")
    println(io, "├─ total samples = $(nsample(samp))")
    println(io, "├─ cores per sample = $(ncore(samp))")
    println(io, "├─ total cores = $(ncore(samp)*nsample(samp))")
    println(io, "└─ bounding box")
    println(
        io,
        "  ├─ min = $(points |> vec |> PointSet |> boundingbox |> minimum)",
    )
    println(
        io,
        "  └─ max = $(points |> vec |> PointSet |> boundingbox |> maximum)",
    )
end

"""$(TYPEDSIGNATURES)
Executes a sample plan, filling an existing `CoreSet` with the results."""
function executeplan!(
    samp::Matrix{Point{2,𝒯}},
    plan::SamplePlan{𝒯};
    stencil::AbstractStencil{𝒯}=SingleCoreStencil{𝒯}(),
    plannedjitter::AbstractJitter{𝒯}=NoJitter{𝒯}(),
    samplerjitter::AbstractJitter{𝒯}=NoJitter{𝒯}(),
    corejitter::AbstractJitter{𝒯}=NoJitter{𝒯}(),
)::Nothing where {𝒯}
    @assert size(samp, 2) == nsample(plan)
    @assert size(samp, 1) == ncore(stencil) "the stencil must have the same number of cores as the sampling matrix has rows, but found $(ncore(stencil)) and $(size(samp,1)), respectively"
    @unpack location, points = plan
    # place points based on the location jitter
    r = location |> first # starting with the first location index and the indices are pre-sorted
    p = points |> first |> plannedjitter
    for j ∈ eachindex(points)
        # if/when the location index changes, update the target point for the core group (each column of the samp matrix)
        if location[j] != r
            r = location[j]
            p = plannedjitter(points[j])
        end
        # apply sampler jitter for the target location (all cores)
        samp[:, j] .= samplerjitter(p)
        # break cores into composite cores (with the stencil) and apply core jitter
        for i ∈ 1:ncore(stencil)
            samp[i, j] = stencil(samp[i, j], i) |> corejitter
        end
    end
    nothing
end

function executeplan!(samp::CoreSet, plan::SamplePlan; kwargs...)::Nothing
    executeplan!(samp.points, plan; kwargs...)
    nothing
end

"""$(TYPEDSIGNATURES)
Executes as sample plan and returns a `CoreSet`"""
function executeplan(
    plan::SamplePlan{𝒯};
    stencil::AbstractStencil{𝒯}=SingleCoreStencil{𝒯}(),
    kwargs...,
)::CoreSet{𝒯} where {𝒯}
    # allocate and pass to the in-place function
    samp = CoreSet(plan, ncore(stencil))
    executeplan!(samp, plan; stencil=stencil, kwargs...)
    return samp
end

"""$(TYPEDSIGNATURES)

Creates a tabular representation of the core set using the [`GeoTable`](https://juliaearth.github.io/GeoStatsDocs/stable/data.html) format. The locations of each *core* in the `CoreSet` are tabulated along with the metadata for each group from the sample plan."""
function GeoTable(samp::CoreSet, plan::SamplePlan)
    @unpack location, round, time, control, points = plan
    C = ncore(samp)
    cols = (
        location=repeat(location, inner=C),
        round=repeat(round, inner=C),
        time=repeat(time, inner=C),
        control=repeat(control, inner=C),
        core=repeat(1:C, nsample(samp)) .|> UInt16,
    )
    georef(cols, samp.points |> vec)
end

end
