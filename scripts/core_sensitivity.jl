using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Monty
using Random
using Distributions
using LaTeXStrings
using LinearAlgebra: Symmetric
using ForwardDiff: gradient
using GlobalSensitivity
using QuasiMonteCarlo

##

import CairoMakie as Mke

##

const plotdir = joinpath(@__DIR__, "plots")
if !isdir(plotdir)
    mkpath(plotdir)
end

##

function coreconcentration(x)
    # the cross-sectional area is always 1, doesn't matter
    core = mixing(x[1], x[2], 1.0, x[3:end]...)
    return core[:analyte]
end

##

γ = 0.9
d = 0.1
Q = 3.0
ρf = 3e3
cf = 0.05
ρs = 1e3
cs = 0.003
𝓁 = 0.5
ℒ = 0.5

x = [γ, d, Q, ρf, cf, ρs, cs, 𝓁, ℒ]

labels = [
    L"$\gamma$",
    L"$d$",
    L"$Q$",
    L"$ρ_f$",
    L"$c_f$",
    L"$ρ_s$",
    L"$c_s$",
    L"$l$",
    L"$L$",
]

g = gradient(coreconcentration, x)

r = 1e6 * 0.01x .* g

fig = Mke.Figure(size=(400, 300))

ax = Mke.Axis(
    fig[1, 1],
    title="Core Concentration Sensitivity",
    xlabel=L"$\frac{u}{100} \cdot \frac{\partial c}{\partial u}$ [ppm]",
    ylabel=L"model parameter ($u$)",
    xlabelsize=18,
    ylabelsize=18,
    yticks=(1:length(r), labels),
)

Mke.barplot!(
    ax,
    1:length(r) |> collect,
    r,
    color=:grey,
    direction=:x,
    strokecolor=:black,
    strokewidth=1,
    width=0.9,
)
Mke.vlines!(ax, 0, color=:black)
Mke.save(joinpath(plotdir, "gradient_sensitivity.png"), fig)
fig

##

function constantsoil(x)
    # the cross-sectional area is always 1, doesn't matter
    core = mixing(x[1], x[2], 1.0, x[3], x[4], x[5], x[6], 0.003, x[7], x[8])
    return core[:analyte]
end

##

samples = 1_000_000
lb = [0.25, 0.05, 0.0, 1e3, 0.04, 500, 0, 0]
ub = [1, 0.2, 5.6, 3e3, 0.09, 1500, 1, 1]

labels =
    [L"$\gamma$", L"$d$", L"$Q$", L"$ρ_f$", L"$c_f$", L"$ρ_s$", L"$l$", L"$L$"]

sen = gsa(
    constantsoil,
    Sobol(order=[0, 1, 2]),
    zip(lb, ub) |> collect,
    samples=samples,
)

##

fig = Mke.Figure(size=(400, 500))

ax = Mke.Axis(
    fig[2, 1],
    xticks=(1:length(labels), labels),
    yticks=xticks = (1:length(labels), labels),
    title="Second Order Sobol Indices",
)
hm = Mke.heatmap!(
    ax,
    1:length(labels),
    1:length(labels),
    (sen.S2 |> Symmetric),
    colormap=:Greys,
    nan_color=:gray,
    colorrange=(0.0, sen.S2 |> maximum),
)
Mke.Colorbar(fig[2, 2], hm)

ax = Mke.Axis(
    fig[1, :],
    xticks=(1:length(labels), labels),
    title="First Order Sobol Indices",
)
Mke.barplot!(ax, 1:length(labels), sen.ST, color=:gray)
Mke.hlines!(ax, 0, color=:black)

Mke.save(joinpath(plotdir, "sobol_sensitivity.png"), fig)
fig

##
