module Gaussian

using UnPack
using LinearAlgebra: Cholesky, Symmetric, LowerTriangular, cholesky!, mul!
using Random: AbstractRNG, randn!
using Distributions
using Meshes: Point
using GeoStatsFunctions: Covariance, Variogram, pairwise, pairwise!

import Base: show
import Random: rand, rand!
import Distributions: MvNormal

using ..Sampling

import ..analytes, ..nanalyte, ..ncore, ..nsample

export updategaussian!, rand!, rand

using DocStringExtensions

#------------------------------------------------------------------------------
export MvNormal

function MvNormal(
    μ::Vector{𝒯},
    σ::Covariance,
    points::AbstractArray{Point{2,𝒯}},
) where {𝒯}
    MvNormal(μ, pairwise(σ, points, points) .|> 𝒯)
end

function MvNormal(μ::Vector{𝒯}, σ::Covariance, samp::CoreSet{𝒯}) where {𝒯}
    MvNormal(μ, σ, samp.points)
end

function MvNormal(μ::𝒯, σ::Covariance, samp::CoreSet{𝒯}) where {𝒯}
    MvNormal(fill(μ, ncore(samp) * nsample(samp)), σ, samp)
end

#------------------------------------------------------------------------------
export GaussianSimulator

"""The `GaussianSimulator` type is for simulating a multivariate distribution over points in space and efficiently reloading the correlation matrix with a new set of points. See the [Efficient Simulation](efficiency.md) example.

## Constructors
    GaussianSimulator(N::Integer, μ::AbstractVector{𝒯}) where {𝒯}

    GaussianSimulator(samp::CoreSet{𝒯}, μ::AbstractVector{𝒯}) where {𝒯}

    GaussianSimulator(samp::CoreSet{𝒯}, μ::𝒯) where {𝒯}
"""
struct GaussianSimulator{𝒯}
    N::UInt16
    μ::Vector{𝒯}
    chol::Matrix{𝒯}
    z::Vector{𝒯}
end

function GaussianSimulator(N::Integer, μ::AbstractVector{𝒯}) where {𝒯}
    GaussianSimulator(
        N |> UInt16,
        μ |> collect,
        Matrix{𝒯}(undef, N, N),
        Vector{𝒯}(undef, N),
    )
end

function GaussianSimulator(samp::CoreSet{𝒯}, μ::AbstractVector{𝒯}) where {𝒯}
    GaussianSimulator(ncore(samp) * nsample(samp), μ)
end

function GaussianSimulator(samp::CoreSet{𝒯}, μ::𝒯) where {𝒯}
    N = ncore(samp) * nsample(samp)
    GaussianSimulator(N, fill(μ, N))
end

"""$(TYPEDSIGNATURES)
Updates the covariance matrix for the `Simulator`` using a `Covariance` model and a group of points."""
function updategaussian!(
    gs::GaussianSimulator{𝒯},
    points::AbstractArray{Point{2,𝒯}},
    c::𝒞,
)::Nothing where {𝒯,𝒞<:Covariance}
    @unpack N, chol = gs
    @assert N == length(points)
    for i ∈ 1:N
        for j ∈ 1:i-1
            chol[i, j] = c(points[i], points[j])
            chol[j, i] = chol[i, j]
        end
        chol[i, i] = c(points[i], points[i])
    end
    cholesky!(chol)
    nothing
end

function GaussianSimulator(
    points::AbstractArray{Point{2,𝒯}},
    μ::AbstractVector{𝒯},
    c::Covariance,
) where {𝒯}
    @assert length(points) == length(μ)
    gs = GaussianSimulator(length(points), μ)
    updategaussian!(gs, points, c)
    return gs
end

function GaussianSimulator(
    points::AbstractArray{Point{2,𝒯}},
    μ::𝒯,
    c::Covariance,
) where {𝒯}
    GaussianSimulator(points, fill(μ, length(points)), c)
end

"""$(TYPEDSIGNATURES)
Samples from the `Simulator` in-place"""
function rand!(
    rng::AbstractRNG,
    gs::GaussianSimulator{𝒯},
    y::AbstractVector{𝒯},
)::Nothing where {𝒯}
    @unpack N, μ, chol, z = gs
    @assert length(y) == N
    randn!(rng, z)
    mul!(y, chol |> Symmetric |> LowerTriangular, z)
    y .+= μ
    nothing
end

"""$(TYPEDSIGNATURES)
Samples from the `Simulator`, returning a new vector"""
function rand(rng::AbstractRNG, gs::GaussianSimulator{𝒯}) where {𝒯}
    y = Vector{𝒯}(undef, gs.N)
    rand!(rng, gs, y)
    return y
end

#------------------------------------------------------------------------------
export GaussianCosimulator

"""The `GaussianCosimulator` type is for simulating a multivariate distribution over points for *two fields simultaneously*, with cross-correlation, and efficiently reloading the correlation matrix with a new set of points. See the [Efficient Simulation](efficiency.md) example.

## Constructors

    GaussianCosimulator(
        N::Integer,
        μ::NamedTuple{𝒦,NTuple{2,Vector{𝒯}}},
        ρ::𝒯,
    )

    GaussianCosimulator(
        samp::CoreSet{𝒯},
        μ::NamedTuple{𝒦,NTuple{2,Vector{𝒯}}},
        ρ::𝒯,
    )

    GaussianCosimulator(
        samp::CoreSet{𝒯},
        μ::NamedTuple{𝒦,NTuple{2,𝒯}},
        ρ::𝒯,
    )
"""
struct GaussianCosimulator{𝒯,𝒦}
    N::UInt16
    ρ::𝒯
    μ::NamedTuple{𝒦,NTuple{2,Vector{𝒯}}}
    σ::Matrix{𝒯}
    chol::Matrix{𝒯}
    z::Vector{𝒯}
    t::Vector{𝒯}
end

function GaussianCosimulator(
    N::Integer,
    μ::NamedTuple{𝒦,NTuple{2,Vector{𝒯}}},
    ρ::𝒯,
) where {𝒯,𝒦}
    @assert zero(𝒯) <= abs(ρ) < one(𝒯) "the correlation coefficient between fields must be between 0 and 1, but got $ρ"
    GaussianCosimulator(
        N |> UInt16,
        ρ,
        μ,
        Matrix{𝒯}(undef, N, N),
        Matrix{𝒯}(undef, 2N, 2N),
        Vector{𝒯}(undef, 2N),
        Vector{𝒯}(undef, 2N),
    )
end

function GaussianCosimulator(
    samp::CoreSet{𝒯},
    μ::NamedTuple{𝒦,NTuple{2,Vector{𝒯}}},
    ρ::𝒯,
) where {𝒯,𝒦}
    GaussianCosimulator(ncore(samp) * nsample(samp), μ, ρ)
end

function GaussianCosimulator(
    samp::CoreSet{𝒯},
    μ::NamedTuple{𝒦,NTuple{2,𝒯}},
    ρ::𝒯,
) where {𝒯,𝒦}
    N = ncore(samp) * nsample(samp)
    GaussianCosimulator(N, ntuple(i -> fill(μ[i], N), 2) |> NamedTuple{𝒦}, ρ)
end

function covariancematrix!(
    σ::Matrix{𝒯},
    points::AbstractArray{Point{2,𝒯}},
    c::Covariance,
)::Nothing where {𝒯}
    @assert size(σ, 1) == size(σ, 2) == length(points)
    N = length(points)
    for i ∈ 1:N
        for j ∈ 1:i-1
            σ[i, j] = c(points[i], points[j])
            σ[j, i] = σ[i, j]
        end
        σ[i, i] = c(points[i], points[i])
    end
    nothing
end

function choleskymatrix!(chol::Matrix{𝒯}, σ::Matrix{𝒯}, ρ::𝒯)::Nothing where {𝒯}
    @assert -1 <= ρ <= 1
    @assert 2 * size(σ, 1) == size(chol, 1)
    @assert 2 * size(σ, 2) == size(chol, 2)
    N = size(σ, 1)
    view(chol, 1:N, 1:N) .= σ
    view(chol, N+1:2N, N+1:2N) .= σ
    view(chol, 1:N, N+1:2N) .= σ
    view(chol, 1:N, N+1:2N) .*= ρ
    view(chol, N+1:2N, 1:N) .= σ
    view(chol, N+1:2N, 1:N) .*= ρ
    cholesky!(chol)
    nothing
end

"""$(TYPEDSIGNATURES)
Updates the covariance matrix for the `Cosimulator` using a `Covariance` model and a group of points."""
function updategaussian!(
    gc::GaussianCosimulator{𝒯},
    points::AbstractArray{Point{2,𝒯}},
    c::Covariance,
) where {𝒯}
    @unpack ρ, σ, chol = gc
    covariancematrix!(σ, points, c)
    choleskymatrix!(chol, σ, ρ)
end

function GaussianCosimulator(
    points::AbstractArray{Point{2,𝒯}},
    μ::NamedTuple{𝒦,NTuple{2,Vector{𝒯}}},
    c::Covariance,
    ρ::𝒯,
) where {𝒯,𝒦}
    @assert zero(𝒯) <= abs(ρ) < one(𝒯) "the correlation coefficient between fields must be between 0 and 1, but got $ρ"
    @assert length(points) == length(μ[1]) == length(μ[2]) "all μ vectors must have the same length and have the same length as the points vector"
    N = length(points)
    gc = GaussianCosimulator(N, μ, ρ)
    updategaussian!(gc, points, c)
    return gc
end

function GaussianCosimulator(
    points::AbstractArray{Point{2,𝒯}},
    μ::NamedTuple{𝒦,NTuple{2,𝒯}},
    c::Covariance,
    ρ::𝒯,
) where {𝒯<:Number,𝒦}
    n = length(points)
    GaussianCosimulator(
        points,
        ntuple(i -> fill(μ[i], n), 2) |> NamedTuple{𝒦},
        c,
        ρ,
    )
end

function GaussianCosimulator(
    samp::CoreSet{𝒯},
    μ::NamedTuple{𝒦,NTuple{2,𝒯}},
    c::Covariance,
    ρ::𝒯,
) where {𝒯<:Number,𝒦}
    GaussianCosimulator(samp.points, μ, c, ρ)
end

analytes(::GaussianCosimulator{𝒯,𝒦}) where {𝒯,𝒦} = 𝒦

nanalyte(::GaussianCosimulator) = 2

ncore(gc::GaussianCosimulator) = gc.N

function show(io::IO, gc::GaussianCosimulator{𝒯,𝒦}) where {𝒯,𝒦}
    println(io, "GaussianCosimulator on analytes $𝒦")
    println(io, "├─ $(ncore(gc)) locations")
    println(io, "└─ cross correlation coefficient = $(gc.ρ)")
end

"""$(TYPEDSIGNATURES)
Samples from the `Cosimulator` in-place, where `y` is a named tuple that includes the two vectors to fill. The keys of the named tuple must include the analytes that the `Cosimulator` expects."""
function rand!(
    rng::AbstractRNG,
    gc::GaussianCosimulator{𝒯,𝒦},
    y::NamedTuple{𝒦},
)::Nothing where {𝒯,𝒦}
    @assert length(y[1]) == length(y[2]) == ncore(gc)
    @unpack N, μ, chol, z, t = gc
    #standard normal starting vector
    randn!(rng, z)
    #apply cholesky L matrix
    mul!(t, chol |> Symmetric |> LowerTriangular, z)
    #copy into the output arrays
    copyto!(y[1], view(t, 1:N))
    copyto!(y[2], view(t, N+1:2N))
    #add mean values
    for i ∈ 1:N
        y[1][i] += μ[1][i]
        y[2][i] += μ[2][i]
    end
    nothing
end

"""$(TYPEDSIGNATURES)
Samples from the `Cosimulator`, returning a new vector"""
function rand(
    rng::AbstractRNG,
    gc::GaussianCosimulator{𝒯,𝒦},
)::NamedTuple{𝒦,NTuple{2,Vector{𝒯}}} where {𝒯,𝒦}
    y = (zeros(𝒯, ncore(gc)), zeros(𝒯, ncore(gc))) |> NamedTuple{𝒦}
    rand!(rng, gc, y)
    return y
end

end
