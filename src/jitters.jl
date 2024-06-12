module Jitters

using Random: AbstractRNG, Xoshiro, rand
using Distributions: UnivariateDistribution, Normal, Beta, mean
using Meshes:
    Point, CartesianGrid, Geometry, GeometrySet, spacing, vertices, coordinates
using StatsAPI: params
using UnPack: @unpack

using DocStringExtensions

#------------------------------------------------------------------------------
export jitterpoint, jitterpoints!

"""$(TYPEDSIGNATURES)
Applies random, isotropic normally distributed jitter to a point in two-dimensional space, where `σ` is the absolute standard deviation of the jitter."""
function jitterpoint(rng::AbstractRNG, point::Point{𝒩,𝒯}, σ::𝒯) where {𝒩,𝒯}
    ntuple(i -> getindex(coordinates(point), i) + σ * randn(rng), 𝒩) |> Point
end

"""$(TYPEDSIGNATURES)
Applies random, isotropic normally distributed jitter to a vector of points, where `σ` is the absolute standard deviation of the jitter."""
function jitterpoints!(
    rng::AbstractRNG,
    points::AbstractVector{Point{2,𝒯}},
    σ::𝒯,
)::Nothing where {𝒯}
    for i ∈ eachindex(points)
        points[i] = jitterpoint(rng, points[i], σ)
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Applies random, isotropic normally distributed jitter to a vector of points *in place*, filling the `jittered` vector with points from `original` after jittering them, where `σ` is the absolute standard deviation of the jitter."""
function jitterpoints!(
    rng::AbstractRNG,
    jittered::AbstractVector{Point{2,𝒯}},
    original::AbstractVector{Point{2,𝒯}},
    σ::𝒯,
)::Nothing where {𝒯}
    @assert length(jittered) == length(original)
    for i ∈ eachindex(original)
        jittered[i] = jitterpoint(rng, original[i], σ)
    end
    nothing
end

#--------------------------------------
export jittergrid

"""$(TYPEDSIGNATURES)
Generates a vector of points where each point is assigned to a random location in each cell of a `CartesianGrid`. The `centering` parameter controls how closely the points tend to be to each cell's centroid. For uniform randomness (no centering) use `centering=1.0`. Higher values induce stronger centering."""
function jittergrid(
    rng::AbstractRNG,
    grid::CartesianGrid{2,𝒯},
    centering::Real,
)::Vector{Point{2,𝒯}} where {𝒯}
    Δx, Δy = spacing(grid)
    β = Beta(centering |> 𝒯, centering |> 𝒯)
    map(grid) do quadrangle
        x₀, y₀ = quadrangle |> vertices |> first |> coordinates
        x = x₀ + Δx * rand(rng, β)
        y = y₀ + Δy * rand(rng, β)
        Point(x, y)
    end
end

#--------------------------------------
export AbstractJitter

abstract type AbstractJitter{𝒯} end

#--------------------------------------
export NoJitter

"""An empty type representing no jitter at all

    NoJitter()
"""
struct NoJitter{𝒯} <: AbstractJitter{𝒯} end

NoJitter(T::Type=Float64) = NoJitter{T}()

(N::NoJitter)(point) = point

#--------------------------------------
export Jitter

"""The `Jitter` type applies random isotropic noise to a `Point`.

    Jitter(d::UnivariateDistribution; seed=Union{Nothing,Integer}=nothing)

The first constructor above accepts a pre-specified distribution, which represents the noise applied in each dimension.

    Jitter(σ; seed::Union{Nothing,Integer}=nothing)

The second constructor assumes a normal distribution for the noise, and receives the standard deviation of that distribution. In both cases a `seed` can be given, which is used to draw random values."""
struct Jitter{𝒯,𝒟<:UnivariateDistribution} <: AbstractJitter{𝒯}
    rng::Xoshiro
    d::𝒟
end

function Jitter(
    d::𝒟;
    seed::Union{Nothing,Integer}=nothing,
) where {𝒟<:UnivariateDistribution}
    @assert d |> mean |> iszero "Jitter distribution must have zero mean, but got $(mean(d))"
    𝒯 = d |> eltype
    Jitter{𝒯,𝒟}(Xoshiro(seed), d)
end

function Jitter(σ::𝒯; seed::Union{Nothing,Integer}=nothing) where {𝒯<:Real}
    @assert σ >= zero(𝒯)
    Jitter(Normal(zero(𝒯), σ), seed=seed)
end

(J::Jitter)(x::Number) = x + rand(J.rng, J.d)

function (J::Jitter{𝒯})(point::Point{𝒩,𝒯}) where {𝒩,𝒯}
    ntuple(i -> J(getindex(coordinates(point), i)), 𝒩) |> Point{𝒩,𝒯}
end

function show(io::IO, J::Jitter{𝒯}) where {𝒯}
    println(io, "Jitter{$𝒯}\n└─ distribution = $(J.d)")
end

#--------------------------------------
export GridCentroidJitter

"""Random variation of points within the dimensions of grid cells. The """
struct GridCentroidJitter{𝒯} <: AbstractJitter{𝒯}
    rng::Xoshiro
    β::Beta{𝒯}
    Δx::𝒯
    Δy::𝒯
end

"""$(TYPEDSIGNATURES)
"""
function GridCentroidJitter(
    spacing::𝒯,
    centering::𝒯;
    seed::Union{Nothing,Integer}=nothing,
) where {𝒯}
    GridCentroidJitter(
        Xoshiro(seed),
        Beta(centering, centering),
        spacing,
        spacing,
    )
end

"""$(TYPEDSIGNATURES)
"""
function GridCentroidJitter(
    Δx::𝒯,
    Δy::𝒯,
    centering::𝒯;
    seed::Union{Nothing,Integer}=nothing,
) where {𝒯}
    GridCentroidJitter(Xoshiro(seed), Beta(centering, centering), Δx, Δy)
end

"""$(TYPEDSIGNATURES)
"""
function GridCentroidJitter(
    grid::CartesianGrid{2,𝒯},
    centering::𝒯;
    seed::Union{Nothing,Integer}=nothing,
) where {𝒯}
    GridCentroidJitter(grid.spacing[1], grid.spacing[2], centering, seed=seed)
end

"""$(TYPEDSIGNATURES)
"""
function GridCentroidJitter(
    grid::CartesianGrid{2,𝒯};
    seed::Union{Nothing,Integer}=nothing,
) where {𝒯}
    GridCentroidJitter(grid, one(𝒯), seed=seed)
end

function (J::GridCentroidJitter{𝒯})(point::Point{2,𝒯})::Point{2,𝒯} where {𝒯}
    @unpack rng, β, Δx, Δy = J
    c = coordinates(point)
    x = c[1] + Δx * rand(rng, β) - Δx / 2
    y = c[2] + Δy * rand(rng, β) - Δy / 2
    Point(x, y)
end

function show(io::IO, J::GridCentroidJitter{𝒯}) where {𝒯}
    @unpack Δx, Δy, β = J
    c = params(β) |> first
    println(
        io,
        "GridCentroidJitter{$𝒯}\n├─ Δx = $Δx\n├─ Δy = $Δy\n└─ centering = $c",
    )
end

end
