module Simulating

using Unitful: @u_str
using Meshes: Point, coordinates
using GeoTables: GeoTable
using GeoStatsFunctions: Covariance, pairwise
using Random: AbstractRNG, Xoshiro, rand, rand!
using Statistics: mean, std
using UnPack: @unpack
using GeoTables: georef, domain
using Printf
using Crayons: Crayon, @crayon_str
using DimensionalData
using Dates: now
using NCDatasets
using CSV
using DataFrames: DataFrame, unstack
using ProgressMeter: Progress, next!
using Distributions:
    Distribution,
    UnivariateDistribution,
    Normal,
    Exponential,
    Uniform,
    TriangularDist

import Base: show
import Statistics: mean
import GeoTables: GeoTable

using ..Checks
using ..MixingModel
using ..Sampling
using ..LeachingModels
using ..Gaussian
using ..Measuring

import ..ncore, ..nsample, ..nanalyte, ..analytes

using DocStringExtensions

#------------------------------------------------------------------------------
export moisturefraction, moistureratio

"""$(TYPEDSIGNATURES)
Computes the moisture fraction (``\\Phi``) from the moisture ratio (``\\Omega``). The "ratio" is the water mass divided by the mass of dry material. The "fraction" is the water mass divided by the total mass. Evaluates
```math
\\Omega / (1 + \\Omega)
```"""
function moisturefraction(Ω)
    @assert isnonnegative(Ω, u"kg/kg") "the moisture ratio cannot be negative, but got $Ω"
    Ω / (1 + Ω)
end

"""$(TYPEDSIGNATURES)
Computes the moisture ratio (``\\Omega``) from the moisture fraction (``\\Phi``). The "ratio" is the water mass divided by the mass of dry material. The "fraction" is the water mass divided by the total mass. Evaluates
```math
\\Phi / (1 - \\Phi)
```"""
function moistureratio(Φ)
    @assert isfractional(Φ, u"kg/kg") "the moisture fraction must be fractional but got $Φ"
    @assert islessthanone(Φ, u"kg/kg") "the moisture ratio should not be computed from a moisture fraction equal to 1 because the result is infinite"
    Φ / (1 - Φ)
end

#------------------------------------------------------------------------------
export Simulation, clear!

"""The `Simulation` type is a container for running simulations. See the [`simulationstack`](@ref) function for more detail.

The fields are
$(TYPEDFIELDS)
"""
struct Simulation{𝒯,𝒩,𝒦}
    # simulation size
    nsamp::UInt16
    ncore::UInt16
    # mixing model inputs
    γ::Matrix{𝒯}
    d::Matrix{𝒯}
    a::𝒯
    Q::Matrix{𝒯}
    ρf::Matrix{𝒯}
    cf::NamedTuple{𝒦,NTuple{𝒩,Matrix{𝒯}}}
    ρs::Matrix{𝒯}
    cs::NamedTuple{𝒦,NTuple{𝒩,Matrix{𝒯}}}
    𝓁::NamedTuple{𝒦,NTuple{𝒩,Matrix{𝒯}}}
    ℒ::Matrix{𝒯}
    # synthetic data at various stages
    cores::Matrix{Sample{𝒩,𝒦,𝒯,𝒯}}
    composites::Vector{Sample{𝒩,𝒦,𝒯,𝒯}}
    measurements::Vector{Measurement{𝒩,𝒦,𝒯}}
end

function Simulation(
    analytes::NTuple{𝒩,Symbol},
    nsamp::Integer,
    ncore::Integer,
    T::Type=Float64;
    a=3e-4, # circle with about 1 cm radius
) where {𝒩}
    Simulation(
        nsamp |> UInt16,
        ncore |> UInt16,
        fill(T(NaN), ncore, nsamp),
        fill(T(NaN), ncore, nsamp),
        convert(T, a),
        fill(T(NaN), ncore, nsamp),
        fill(T(NaN), ncore, nsamp),
        NamedTuple{analytes}(ntuple(_ -> fill(T(NaN), ncore, nsamp), 𝒩)),
        fill(T(NaN), ncore, nsamp),
        NamedTuple{analytes}(ntuple(_ -> fill(T(NaN), ncore, nsamp), 𝒩)),
        NamedTuple{analytes}(ntuple(_ -> fill(T(NaN), ncore, nsamp), 𝒩)),
        fill(T(NaN), ncore, nsamp),
        Matrix{Sample{𝒩,analytes,T,T}}(undef, ncore, nsamp),
        Vector{Sample{𝒩,analytes,T,T}}(undef, nsamp),
        Vector{Measurement{𝒩,analytes,T}}(undef, nsamp),
    )
end

function Simulation(
    analytes::NTuple{𝒩,Symbol},
    samp::CoreSet{𝒯};
    a=1.25e-3,
) where {𝒩,𝒯}
    Simulation(analytes, nsample(samp), ncore(samp), 𝒯, a=(a |> 𝒯))
end

nsample(S::Simulation) = S.nsamp

ncore(S::Simulation) = S.ncore

nanalyte(::Simulation{𝒯,𝒩}) where {𝒯,𝒩} = 𝒩

analytes(::Simulation{𝒯,𝒩,𝒦}) where {𝒯,𝒩,𝒦} = 𝒦

analytes(::Simulation{𝒯,𝒩,𝒦}, i::Integer) where {𝒯,𝒩,𝒦} = 𝒦[i]

"""$(TYPEDSIGNATURES)
Clears all fields of a `Simulation`"""
function clear!(sim::Simulation{𝒯,𝒩})::Nothing where {𝒯,𝒩}
    sim.γ .= NaN |> 𝒯
    sim.d .= NaN |> 𝒯
    sim.Q .= NaN |> 𝒯
    sim.ρf .= NaN |> 𝒯
    sim.ρs .= NaN |> 𝒯
    sim.ℒ .= NaN |> 𝒯
    for i ∈ 1:nanalyte(sim)
        sim.cf[i] .= NaN |> 𝒯
        sim.cs[i] .= NaN |> 𝒯
        sim.𝓁[i] .= NaN |> 𝒯
    end
    for i ∈ eachindex(sim.cores)
        sim.cores[i] = zero(sim.cores[1])
    end
    for i ∈ eachindex(sim.composites)
        sim.composites[i] = zero(sim.composites[1])
    end
    for i ∈ eachindex(sim.measurements)
        sim.measurements[i] = zero(sim.measurements[1])
    end
    nothing
end

function description(x)
    μ = @sprintf "%.3g" mean(x)
    σ = @sprintf "%.3g" std(x)
    rσ = @sprintf "%.1f" (100 * std(x) / abs(mean(x)))
    lo = @sprintf "%.3g" minimum(x)
    locrayon = if (minimum(x) < 0)
        Crayon(foreground=:red, bold=true, italics=true)
    else
        Crayon(foreground=:cyan)
    end
    hi = @sprintf "%.3g" maximum(x)
    if isnan(mean(x))
        "$μ ± $σ ($rσ %) ∈ [$lo, $hi]" |> crayon"yellow"
    else
        d = string("$μ" |> Crayon(bold=true))
        d *= " ± "
        d *= string("$σ" |> Crayon(foreground=:dark_gray, italics=true))
        d *= " ($rσ %)" |> Crayon(foreground=:light_blue, bold=true) |> string
        d *= " ∈ ["
        d *= string("$lo" |> locrayon)
        d *= ", "
        d *= string("$hi" |> Crayon(foreground=:magenta))
        d *= "]"
        return d
    end
end

function show(io::IO, sim::Simulation{𝒯,𝒩,𝒦}) where {𝒯,𝒩,𝒦}
    s = 𝒩 > 1 ? "analytes" : "analyte"
    k = join(map(string, 𝒦), ", ")
    println(io, "Simulation{$𝒯} with $𝒩 $s: $k")
    println(io, "├─ total samples = $(nsample(sim))")
    println(io, "├─ cores per sample = $(ncore(sim))")
    println(io, "├─ total cores = $(nsample(sim) * ncore(sim))")
    println(io, "└─ simulated parameters μ ± σ (rsd) ∈ [max, min]")
    println(io, "  ├─ core cross-sectional area (a)\n  │    $(sim.a)")
    println(io, "  ├─ applied feedstock (Q)\n  │    ", description(sim.Q))
    println(io, "  ├─ feedstock fraction (γ)\n  │    ", description(sim.γ))
    println(io, "  ├─ core depth (d)\n  │    ", description(sim.d))
    println(
        io,
        "  ├─ feedstock bulk density (ρf)\n  │    ",
        description(sim.ρf),
    )
    println(io, "  ├─ feedstock concentrations (cf)")
    for i ∈ 1:𝒩
        println(io, "  │    $(𝒦[i]): ", description(sim.cf[i]))
    end
    println(io, "  ├─ soil bulk density (ρs)\n  │    ", description(sim.ρs))
    println(io, "  ├─ background soil concentrations (cs)")
    for i ∈ 1:𝒩
        println(io, "  │    $(𝒦[i]): ", description(sim.cs[i]))
    end
    println(io, "  ├─ feedstock leached fractions (𝓁)")
    for i ∈ 1:𝒩
        println(io, "  │    $(𝒦[i]): ", description(sim.𝓁[i]))
    end
    println(
        io,
        "  └─ feedstock mass loss fraction (ℒ)\n       ",
        description(sim.ℒ),
    )
end

#------------------------------------------------------------------------------
export spreading!

"""$(TYPEDSIGNATURES)
Samples spreading/application rates from a distribution or `Simulator` (or `Cosimulator`), automatically ignoring control points as identified by the provided `plan`."""
function spreading!(
    rng::AbstractRNG,
    sim::Simulation{𝒯,𝒩,𝒦},
    plan::SamplePlan,
    X::UnivariateDistribution,
)::Nothing where {𝒯,𝒩,𝒦}
    @assert eltype(X) == 𝒯
    @assert isnonnegative(X)
    @assert nsample(sim) == nsample(plan)
    @unpack control = plan
    @unpack Q = sim
    for j ∈ 1:nsample(sim)
        if !control[j]
            rand!(rng, X, view(Q, :, j))
        else
            Q[:, j] .= zero(𝒯)
        end
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Samples any object for which `rand!` can be called to apply feedstock, ignoring controls."""
function spreading!(
    rng::AbstractRNG,
    sim::Simulation{𝒯},
    plan::SamplePlan,
    X,#::MvNormal,
)::Nothing where {𝒯}
    #@assert eltype(X) == 𝒯
    @assert nsample(sim) == nsample(plan)
    @unpack control = plan
    @unpack Q = sim
    rand!(rng, X, view(sim.Q, :))
    for j ∈ 1:nsample(sim)
        if control[j]
            Q[:, j] .= zero(𝒯)
        end
    end
end

#------------------------------------------------------------------------------
export unmixed!, triangularmixing!, uniformmixing!, exponentialmixing!

"""$(TYPEDSIGNATURES)
Sets feedstock fractions and sample depths assuming the feedstock is resting in a layer on top of soil (no mixing)."""
function unmixed!(
    rng::AbstractRNG,
    sim::Simulation{𝒯};
    depth::UnivariateDistribution,
)::Nothing where {𝒯}
    @assert eltype(depth) == 𝒯
    @assert isnonnegative(depth)
    @unpack γ, d = sim
    for i ∈ eachindex(γ)
        d[i] = rand(rng, depth)
        γ[i] = one(𝒯)
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Assumes feedstock mixing profiles with a wedge shape"""
function triangularmixing!(
    rng::AbstractRNG,
    sim::Simulation{𝒯};
    depth::UnivariateDistribution,
    upper::UnivariateDistribution,
)::Nothing where {𝒯}
    @assert eltype(depth) == eltype(upper) == 𝒯
    @assert isnonnegative(depth)
    @assert isnonnegative(upper)
    @unpack γ, d = sim
    for i ∈ eachindex(γ)
        d[i] = rand(rng, depth)
        Γ = TriangularDist(0, rand(rng, upper), 0)
        γ[i] = feedstockfraction(Γ, d[i])
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Assumes feedstock mixing profiles uniform over a depth interval"""
function uniformmixing!(
    rng::AbstractRNG,
    sim::Simulation{𝒯};
    depth::UnivariateDistribution,
    upper::UnivariateDistribution,
)::Nothing where {𝒯}
    @assert eltype(depth) == eltype(upper) == 𝒯
    @assert isnonnegative(depth)
    @assert isnonnegative(upper)
    @unpack γ, d = sim
    for i ∈ eachindex(γ)
        d[i] = rand(rng, depth)
        Γ = Uniform(0, rand(rng, upper))
        γ[i] = feedstockfraction(Γ, d[i])
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Assumes exponentially decaying feedstock mixing profiles"""
function exponentialmixing!(
    rng::AbstractRNG,
    sim::Simulation{𝒯};
    depth::UnivariateDistribution,
    scale::UnivariateDistribution,
)::Nothing where {𝒯}
    #@assert eltype(depth) == eltype(scale) == 𝒯
    @assert isnonnegative(depth)
    @assert isnonnegative(scale)
    @unpack γ, d = sim
    for i ∈ eachindex(γ)
        d[i] = rand(rng, depth)
        Γ = Exponential(rand(rng, scale))
        γ[i] = feedstockfraction(Γ, d[i])
    end
    nothing
end

#------------------------------------------------------------------------------
export feedstockconcentration!

"""$(TYPEDSIGNATURES)
Sets feedstock concentrations to constants"""
function feedstockconcentration!(
    sim::Simulation{𝒯,𝒩,𝒦},
    concentrations::NamedTuple{𝒦,NTuple{𝒩,𝒯}},
)::Nothing where {𝒯,𝒩,𝒦}
    for analyte ∈ 𝒦
        sim.cf[analyte] .= concentrations[k]
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Sets feedstock concentrations using a mean and relative standard deviation for each analyte"""
function feedstockconcentration!(
    rng::AbstractRNG,
    sim::Simulation{𝒯,𝒩,𝒦},
    μ::NamedTuple{𝒦,NTuple{𝒩,𝒯}},
    σ::𝒯,
)::Nothing where {𝒯,𝒩,𝒦}
    for analyte ∈ 𝒦
        X = Normal(μ[analyte], σ * μ[analyte])
        rand!(rng, X, sim.cf[analyte])
    end
    nothing
end

#------------------------------------------------------------------------------
export soilconcentration!

"""$(TYPEDSIGNATURES)
Sets soil analyte concentrations using any distribution"""
function soilconcentration!(
    rng::AbstractRNG,
    sim::Simulation{𝒯},
    analyte::Symbol,
    X::Distribution,
)::Nothing where {𝒯}
    @assert eltype(X) == 𝒯
    rand!(rng, X, view(sim.cs[analyte], :))
    nothing
end

"""$(TYPEDSIGNATURES)
Sets soil analyte concentrationss using means and relative standard deviation"""
function soilconcentration!(
    rng::AbstractRNG,
    sim::Simulation{𝒯,𝒩,𝒦},
    μ::NamedTuple{𝒦,NTuple{𝒩,𝒯}},
    σ::𝒯,
)::Nothing where {𝒯,𝒩,𝒦}
    for analyte ∈ 𝒦
        X = Normal(μ[analyte], σ * μ[analyte])
        soilconcentration!(rng, sim, analyte, X)
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Sets soil analyte concentrations using a multivariate normal distribution defined by a `Covariance`"""
function soilconcentration!(
    rng::AbstractRNG,
    sim::Simulation{𝒯},
    samp::CoreSet{𝒯},
    analyte::Symbol,
    μ,
    σ::Covariance,
)::Nothing where {𝒯}
    @assert ncore(sim) == ncore(samp)
    @assert nsample(sim) == nsample(samp)
    soilconcentration!(rng, sim, analyte, MvNormal(μ, σ, samp))
    nothing
end

"""$(TYPEDSIGNATURES)
Sets soil analyte concentrations using a [`GaussianSimulator`](@ref)"""
function soilconcentration!(
    rng::AbstractRNG,
    sim::Simulation{𝒯},
    gc::GaussianCosimulator{𝒯,𝒦},
)::Nothing where {𝒯,𝒦}
    @unpack cs = sim
    # the simulation may have more analytes than the GC
    k = analytes(gc)
    y = (cs[k[1]], cs[k[2]]) |> NamedTuple{𝒦}
    rand!(rng, gc, y)
    nothing
end

#------------------------------------------------------------------------------
export leaching!

"""$(TYPEDSIGNATURES)
Sets leached fractions for one analyte using a Leaching Model"""
function leaching!(
    𝓁::Matrix{𝒯},
    model::AbstractLeachingModel{𝒯},
    time::AbstractVector{𝒯},
)::Nothing where {𝒯}
    for (i, t) ∈ enumerate(time)
        𝓁[:, i] .= model(t)
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Sets leached fractions for one analyte using a Leaching Model"""
function leaching!(
    sim::Simulation{𝒯,𝒩,𝒦},
    analyte::Symbol,
    model::AbstractLeachingModel{𝒯},
    time::AbstractVector{𝒯},
)::Nothing where {𝒯,𝒩,𝒦}
    @assert nsample(sim) == length(time)
    leaching!(sim.𝓁[analyte], model, time)
    nothing
end

"""$(TYPEDSIGNATURES)
Sets leached fractions using a Leaching Model and times from a [`SamplePlan`](@ref)"""
function leaching!(
    sim::Simulation{𝒯,𝒩,𝒦},
    analyte::Symbol,
    model::AbstractLeachingModel{𝒯},
    plan::SamplePlan{𝒯},
)::Nothing where {𝒯,𝒩,𝒦}
    leaching!(sim, analyte, model, plan.time)
    nothing
end

"""$(TYPEDSIGNATURES)
Sets leached fractions using a Leaching Model for each analyte"""
function leaching!(
    sim::Simulation{𝒯,𝒩,𝒦},
    models::NamedTuple{𝒦,NTuple{𝒩,ℳ}},
    time::AbstractVector{𝒯},
)::Nothing where {𝒯,𝒩,𝒦,ℳ<:AbstractLeachingModel}
    @assert nsample(sim) == length(time)
    for i ∈ 1:𝒩
        leaching!(sim.𝓁[i], models[i], time)
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Sets leached fractions using a Leaching Model for each analyte and times from a [`SamplePlan`](@ref)"""
function leaching!(
    sim::Simulation{𝒯,𝒩,𝒦},
    models::NamedTuple{𝒦,NTuple{𝒩,ℳ}},
    plan::SamplePlan{𝒯},
)::Nothing where {𝒯,𝒩,𝒦,ℳ<:AbstractLeachingModel}
    leaching!(sim, models, plan.time)
    nothing
end

#--------------------------------------
export massloss!

"""$(TYPEDSIGNATURES)
Sets feedstock mass loss fraction using a function of analyte leached fractions"""
function massloss!(
    sim::Simulation{𝒯,𝒩,𝒦},
    plan::SamplePlan{𝒯},
    loss::ℱ,
)::Nothing where {𝒯,𝒩,𝒦,ℱ<:Function}
    @unpack ℒ, 𝓁 = sim
    @unpack time = plan
    for i ∈ 1:ncore(sim), j ∈ 1:nsample(sim)
        if time[j] < 0
            ℒ[i, j] = zero(𝒯)
        else
            𝓁ᵢⱼ = ntuple(k -> 𝓁[k][i, j], 𝒩) |> NamedTuple{𝒦}
            ℒ[i, j] = loss(𝓁ᵢⱼ)
            @assert zero(𝒯) <= ℒ[i, j] <= one(𝒯) "the feedstock mass loss fraction cannot be outside the [0,1] interval but encountered $(ℒ[i,j]) at time $(time[j]) with elemental loss fractions: $𝓁ᵢⱼ"
        end
    end
    nothing
end

#------------------------------------------------------------------------------
export core!

"""$(TYPEDSIGNATURES)
Executes the core collection process internally in a [`Simulation`](@ref)"""
function core!(sim::Simulation{𝒯,𝒩,𝒦})::Nothing where {𝒯,𝒩,𝒦}
    @unpack γ, d, a, Q, ρf, cf, ρs, cs, 𝓁, ℒ, cores = sim
    for j ∈ 1:nsample(sim), i ∈ 1:ncore(sim)
        cores[i, j] = MixingModel.mixing(
            γ[i, j],
            d[i, j],
            a,
            Q[i, j],
            ρf[i, j],
            ntuple(k -> cf[k][i, j], 𝒩) |> NamedTuple{𝒦},
            ρs[i, j],
            ntuple(k -> cs[k][i, j], 𝒩) |> NamedTuple{𝒦},
            ntuple(k -> 𝓁[k][i, j], 𝒩) |> NamedTuple{𝒦},
            ℒ[i, j],
        )
    end
    nothing
end

#------------------------------------------------------------------------------
export composite!

"""$(TYPEDSIGNATURES)
Executes the core compositing process internally in a [`Simulation`](@ref)"""
function composite!(sim::Simulation)::Nothing
    @unpack cores, composites = sim
    for j ∈ 1:nsample(sim)
        composites[j] = sum(view(cores, :, j))
    end
    nothing
end

#------------------------------------------------------------------------------
export measure!

"""$(TYPEDSIGNATURES)
Executes the composite sample measurement process internally in a [`Simulation`](@ref) using a tuple of relative standard deviations for analytes and a relative standard deviation for the mass"""
function measure!(
    rng::AbstractRNG,
    sim::Simulation{𝒯,𝒩,𝒦},
    σ::NamedTuple{𝒦,NTuple{𝒩,𝒯}},
    σₘ::𝒯,
)::Nothing where {𝒯,𝒩,𝒦}
    @unpack composites, measurements = sim
    for i ∈ eachindex(composites)
        measurements[i] = measure(rng, composites[i], σ, σₘ)
    end
    nothing
end

#------------------------------------------------------------------------------
export analyze!

"""$(TYPEDSIGNATURES)
Executes the entire core collection, compositing, and measurement sequence. That is, the `analyze!` function
1. calls [`core!`](@ref)
2. calls [`composite!`](@ref)
3. calls [`measure`](@ref), passing inputs through to that function"""
function analyze!(
    rng::AbstractRNG,
    sim::Simulation{𝒯,𝒩,𝒦},
    σ::NamedTuple{𝒦,NTuple{𝒩,𝒯}},
    args...,
)::Nothing where {𝒯,𝒩,𝒦}
    core!(sim)
    composite!(sim)
    measure!(rng, sim, σ, args...)
    nothing
end

#------------------------------------------------------------------------------

function GeoTable(
    measurements::Vector{Measurement{𝒩,𝒦,𝒯}},
    cols,
    geom,
) where {𝒩,𝒦,𝒯}
    conc = ntuple(𝒩) do i
        getindex.(measurements, i)
    end |> NamedTuple{𝒦}
    mass = getfield.(measurements, :mass)
    georef(merge(cols, conc, (mass=mass,)), geom)
end

"""$(TYPEDSIGNATURES)

Creates a tabular representation of the simulation using the [`GeoTable`](https://juliaearth.github.io/GeoStatsDocs/stable/data.html) format. The measurement for each composited location is assigned from the sample plan."""
function GeoTable(sim::Simulation, plan::SamplePlan)
    @assert nsample(sim) == nsample(plan)
    @unpack measurements = sim
    gt = Sampling.GeoTable(plan)
    cols = values(gt)
    geom = domain(gt)
    GeoTable(measurements, cols, geom)
end

function mean(x::AbstractVector{Point{𝒩,𝒯}})::Point{𝒩,𝒯} where {𝒩,𝒯}
    sum(p -> coordinates(p) / length(x), x) |> Point
end

"""$(TYPEDSIGNATURES)

Creates a tabular representation of the simulation using the [`GeoTable`](https://juliaearth.github.io/GeoStatsDocs/stable/data.html) format. The measurement for each composited location is assigned by averaging the locations of each core in the core set."""
function GeoTable(sim::Simulation, samp::CoreSet, plan::SamplePlan)
    @assert nsample(sim) == nsample(samp) == nsample(plan)
    @assert ncore(sim) == ncore(samp)
    @unpack measurements = sim
    gt = Sampling.GeoTable(plan)
    cols = values(gt)
    geom = map(mean, eachcol(samp.points))
    GeoTable(measurements, cols, geom)
end

#------------------------------------------------------------------------------
export simulationstack

"""$(TYPEDSIGNATURES)
Runs repeated simulations. See the [Efficient Simulation](efficiency.md) example."""
function simulationstack(
    simulate!::ℱ,
    nrealization::Integer,
    sim::Simulation{𝒯},
    samp::CoreSet{𝒯},
    plan::SamplePlan{𝒯};
    show_progress::Bool=true,
) where {ℱ<:Function,𝒯}
    sims = DimStack((
        data=DimArray(
            Array{𝒯}(undef, nrealization, nanalyte(sim) + 1, nsample(sim)),
            (
                realization=1:nrealization,
                analyte=[analytes(sim)..., (:mass,)...],
                sample=1:nsample(sim),
            ),
        ),
        x=DimArray(
            Matrix{𝒯}(undef, nrealization, nsample(sim)),
            (:realization, :sample),
        ),
        y=DimArray(
            Matrix{𝒯}(undef, nrealization, nsample(sim)),
            (:realization, :sample),
        ),
        control=DimArray(plan.control, (:sample,)),
        location=DimArray(plan.location, (:sample,)),
        round=DimArray(plan.round, (:sample,)),
        time=DimArray(plan.time, (:sample,)),
    ))

    data = sims[:data]
    x = sims[:x]
    y = sims[:y]

    prog = Progress(nrealization, enabled=show_progress)
    for i ∈ 1:nrealization
        # execute the simulation routine, which should fill sim
        simulate!()
        # fill the output arrays with this realization's results
        for (j, analyte) ∈ analytes(sim) |> enumerate
            for k ∈ eachindex(sim.measurements)
                data[i, j, k] = sim.measurements[k][analyte]
            end
        end
        for k ∈ eachindex(sim.measurements)
            data[i, nanalyte(sim)+1, k] = sim.measurements[k].mass
        end
        for j ∈ 1:nsample(sim)
            point = view(samp.points, :, j) |> mean
            x[i, j] = getindex(coordinates(point), 1)
            y[i, j] = getindex(coordinates(point), 2)
        end
        next!(prog)
    end

    return sims
end

#--------------------------------------
export tonetcdf, tocsv

"""$(TYPEDSIGNATURES)
Writes a simulation stack to netcdf format"""
function tonetcdf(fn::AbstractString, sims::DimStack)::Nothing
    T = sims[:data] |> eltype
    ds = NCDataset(fn, "c")
    defDim(ds, "realization", size(sims, 1))
    defDim(ds, "analyte", size(sims, 2))
    defDim(ds, "sample", size(sims, 3))
    # simulated data
    data = defVar(ds, "data", T, ("realization", "analyte", "sample"))
    data[:, :, :] = sims.data
    # average core coordinates
    x = defVar(ds, "x", T, ("realization", "sample"))
    x[:, :] = sims[:x]
    y = defVar(ds, "y", T, ("realization", "sample"))
    y[:, :] = sims[:y]
    # associated information
    control = defVar(ds, "control", Int8, ("sample",))
    control[:] = sims[:control].data .|> Int8
    location = defVar(ds, "location", UInt16, ("sample",))
    location[:] = sims[:location].data
    round = defVar(ds, "round", UInt16, ("sample",))
    round[:] = sims[:round].data
    time = defVar(ds, "time", T, ("sample",))
    time[:] = sims[:time].data
    # attributes
    ds.attrib["analytes"] = string.(sims[:data].dims[2])
    ds.attrib["datetime"] = "file created: $(now() |> string)"
    ds.attrib["comments"] = "simulated geochemical data generated by Monty"
    close(ds)
    nothing
end

"""$(TYPEDSIGNATURES)
Writes one realization of a simulation stack to csv format"""
function tocsv(
    fn::AbstractString,
    sims::DimStack,
    realization::Integer,
)::Nothing
    CSV.write(
        fn,
        unstack(sims[realization, :, :] |> DataFrame, :analyte, :data),
    )
    nothing
end

end
