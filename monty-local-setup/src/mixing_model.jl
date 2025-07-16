module MixingModel

using Distributions: UnivariateDistribution, cdf
using Unitful: @u_str, Quantity, ustrip
using UnPack: @unpack
using Printf

import Base: +, show, length, getindex, zero

using ..Checks

import ..analytes, ..nanalyte

using DocStringExtensions

#------------------------------------------------------------------------------
export feedstockfraction, feedstockmass, feedstockelementmass
export feedstockthickness
export soilmass, soilelementmass, samplemass

"""$(TYPEDSIGNATURES)
Computes the fraction of applied original feedstock contained in a sample, given the mixing profile `Γ`, which must be a univariate distribution, and the sample depth `d`. This function simply wraps a `cdf` evaluation (cumulative probability density) and includes checks for physical consistency."""
function feedstockfraction(Γ::UnivariateDistribution, d)
    @assert isnonnegative(Γ) "feedstock mixing distributions must have strictly non-negative (>= 0) support, but found a minimum support of $(minimum(Γ)) for a $(typeof(Γ)) distribution"
    @assert isnonnegative(d, u"m") "sampling depth must be strictly non-negative (>= 0), but got $d"
    cdf(Γ, d)
end

function feedstockfraction(Γ::UnivariateDistribution, d::Quantity)
    feedstockfraction(Γ, d |> ustrip)
end

"""$(TYPEDSIGNATURES)
Computes the mass of feedstock contained in a sample, given the application rate `Q`, the sampled fraction `γ`, and the feedstock mass loss fraction `ℒ`. Evaluates
```math
Q \\gamma (1 - \\mathscr{L})
```"""
function feedstockmass(Q, γ, ℒ)
    @assert isnonnegative(Q, u"kg/m^2")
    @assert isfractional(γ, u"kg/kg")
    @assert isfractional(ℒ, u"kg/kg")
    Q * γ * (1 - ℒ)
end

"""$(TYPEDSIGNATURES)
Computes the mass of an individual element or analyte in a sample, given the application rate `Q`, the sampled feedstock fraction `γ`, the analyte loss fraction `𝓁`, and the feedstock mass loss fraction `ℒ`. Evaluates
```math
Q \\gamma (1 - \\mathscr{l}) \\mathscr{L}
```"""
function feedstockelementmass(Q, γ, 𝓁, cf)
    @assert isnonnegative(Q, u"kg/m^2")
    @assert isfractional(γ, u"kg/kg")
    @assert isfractional(𝓁, u"kg/kg")
    @assert isfractional(cf, u"kg/kg")
    Q * γ * (1 - 𝓁) * cf
end

"""$(TYPEDSIGNATURES)
Computes the thickness of the feedstock in a sample, given the application rate `Q`, the sampled feedstock fraction `γ`, and the feedstock bulk density in the sample `ρs`. Evaluates
```math
Q \\gamma / \\rho_f
```"""
function feedstockthickness(Q, γ, ρf, ℒ)
    @assert ispositive(ρf, u"kg/m^3")
    feedstockmass(Q, γ, ℒ) / ρf
end

"""$(TYPEDSIGNATURES)
Computes the soil mass in a sample, given the soil bulk density `ρs` and the soil/sample depth `d`. Note that in analytical model, `d` represents an *apparent* soil depth after accounting for the thickness of feedstock material. This function simply performs physical checks and evaluates
```math
\\rho_s d
```"""
function soilmass(ρs, d)
    @assert ispositive(ρs, u"kg/m^3")
    @assert isnonnegative(d, u"m")
    ρs * d
end

"""$(TYPEDSIGNATURES)
Computes the mass of an element/analyte in soil, given the soil bulk density `ρs`, the soil/sample depth `d`, and the concentration of the element in the soil `cs`. Evaluates
```math
\\rho_s d c_s
```"""
function soilelementmass(ρs, d, cs)
    @assert ispositive(ρs, u"kg/m^3")
    @assert isnonnegative(d, u"m")
    @assert isfractional(cs, u"kg/kg") "Soil concentrations must be fractional (between zero and one, inclusive of the boundaries) but got $cs"
    ρs * d * cs
end

"""$(TYPEDSIGNATURES)
Computes the mass of a sample given the sampled soil's mass per unit area `M` and sampled cross-sectional area `a`. Evaluates
```
M a
```"""
function samplemass(M, a)
    @assert isnonnegative(M, u"kg/m^2")
    @assert isnonnegative(a, u"m^2")
    M * a
end

#------------------------------------------------------------------------------
export Sample, singlecore, zero

#=====
Q - 1 per deployment
ρf - function of 𝓁
cf - 1 per deployment, probably just one overall
ρs - 1 overall
cs - 1 overall
d - 1 per deployment, probably just 1 overall
a - 1 per deployment, probably just 1 overall
Γ - 1 per deployment
ℒ - function of 𝓁
𝓁 - 1 for each element and each deployment
=====#

"""$(TYPEDSIGNATURES)
The `Sample` type represents a mass of soil/material with some number of analyte concentrations. When defined, all samples represent an individual soil core, which can be combined by addition (the `+` operator) into composite samples.

The `Sample` type has three fields:

| Field | Description |
| ---: | :--- |
| `concentrations` | a `NamedTuple` containing exact sample concentrations |
| `mass` | the exact mass of the sample |
| `cores` | the number of individual cores combined into the composite sample |

A `Sample` represents the true, exact properties of a piece of soil and is immutable. Once defined, it can't be modifed. When compositing samples (adding them together), the composite concentrations are computed exactly by mass weighted average. Samples are created from soil and deployment parameters using the [`mixing`](@ref) function. To simulate measurement noise/uncertainty, a sample is passed to the `measure` function, which applies gaussian noise defined in terms of relative standard deviations.
"""
struct Sample{𝒩,𝒦,𝒯<:Number,𝒰<:Number}
    concentrations::NamedTuple{𝒦,NTuple{𝒩,𝒯}}
    mass::𝒰
    cores::UInt16
end

# indexing into the sample passes directly through to the concentrations
getindex(sample::Sample, idx) = getindex(sample.concentrations, idx)

"""$(TYPEDSIGNATURES)
Returns the number of analytes in a [`Sample`](@ref)"""
nanalyte(::Sample{𝒩}) where {𝒩} = 𝒩

"""$(TYPEDSIGNATURES)
Returns the names/symbols of the analytes in a [`Sample`](@ref)"""
analytes(::Sample{𝒩,𝒦}) where {𝒩,𝒦} = 𝒦

# mass is extrinsic
# concentrations are intrinsic
# cores are counted
function +(A::Sample{𝒩,𝒦}, B::Sample{𝒩,𝒦}) where {𝒩,𝒦}
    M = A.mass + B.mass
    c = ntuple(𝒩) do j
        (A.mass * A[j] + B.mass * B[j]) / M
    end |> NamedTuple{𝒦}
    cores = A.cores + B.cores
    Sample(c, M, cores)
end

function zero(::Sample{𝒩,𝒦,𝒯,𝒰}) where {𝒩,𝒦,𝒯,𝒰}
    Sample(ntuple(i -> zero(𝒯), 𝒩) |> NamedTuple{𝒦}, zero(𝒰), zero(UInt16))
end

function show(io::IO, sample::Sample{𝒩,𝒦,𝒯}) where {𝒩,𝒦,𝒯}
    @unpack mass, cores = sample
    a = 𝒩 > 1 ? "analytes" : "analyte"
    c = cores > 1 ? "cores" : "core"
    println(io, "Sample{$𝒯} ($𝒩 $a, $cores $c)\nmass = $mass")
    L = map(length ∘ string, 𝒦) |> maximum
    for (n, k) ∈ 𝒦 |> enumerate
        stick = (n < 𝒩) ? "├─" : "└─"
        println(io, "$stick $(rpad(k, L)) $(sample[k])")
    end
end

function show(io::IO, samples::Vector{Sample{𝒩,𝒦,𝒯,𝒰}}) where {𝒩,𝒦,𝒯,𝒰}
    L = length(samples)
    data = zeros(Any, L, 𝒩 + 1)
    for i ∈ 1:𝒩
        data[:, i] = getindex.(samples, 𝒦[i])
    end
    data[:, end] = getfield.(samples, :mass)
    println(io, data)
end

"""$(TYPEDSIGNATURES)
Creates a single-core [`Sample`](@ref) by explicit definition of its concentrations and mass"""
function singlecore(concentrations::NamedTuple, mass)
    Sample(concentrations, mass, one(UInt16))
end

#------------------------------------------------------------------------------
export mixing

"""$(TYPEDSIGNATURES)
Evaluates the mixing and leaching model for a single piece of soil, returning a [`Sample`](@ref)."""
function mixing(
    γ::Number,
    d::Number,
    a::Number,
    Q::Number,
    ρf::Number,
    cf::NamedTuple{𝒦,NTuple{𝒩,𝒯}},
    ρs::Number,
    cs::NamedTuple{𝒦,NTuple{𝒩,𝒯}},
    𝓁::NamedTuple{𝒦,NTuple{𝒩,𝒰}},
    ℒ::Number,
) where {𝒩,𝒦,𝒯,𝒰}
    h = feedstockthickness(Q, γ, ρf, ℒ)
    M = feedstockmass(Q, γ, ℒ) + soilmass(ρs, d - h)
    c = ntuple(𝒩) do j
        cfⱼ = feedstockelementmass(Q, γ, 𝓁[j], cf[j])
        csⱼ = soilelementmass(ρs, d - h, cs[j])
        (cfⱼ + csⱼ) / M
    end |> NamedTuple{𝒦}
    return Sample(c, samplemass(M, a), one(UInt16))
end

"""$(TYPEDSIGNATURES)
Evaluates the mixing and leaching model for a single piece of soil, returning a [`Sample`](@ref)."""
function mixing(
    γ::Number,
    d,
    a,
    Q,
    ρf,
    cf::Number,
    ρs,
    cs::Number,
    𝓁::Number,
    ℒ,
)
    mixing(γ, d, a, Q, ρf, (analyte=cf,), ρs, (analyte=cs,), (analyte=𝓁,), ℒ)
end

function mixing(Γ::UnivariateDistribution, d, args...)
    mixing(feedstockfraction(Γ, d), d, args...)
end

end
