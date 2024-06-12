<div align="center">
<img src="docs/src/assets/logo.png" width="200" height="200">
</div>

# Monty.jl

[![checks](https://github.com/LithosCarbon/Monty.jl/actions/workflows/checks.yml/badge.svg)](https://github.com/LithosCarbon/Monty.jl/actions/workflows/checks.yml) [![benchmarks](https://github.com/LithosCarbon/Monty.jl/actions/workflows/benchmarks.yml/badge.svg?branch=main)](https://github.com/LithosCarbon/Monty.jl/actions/workflows/benchmarks.yml) [![docs](https://github.com/LithosCarbon/Monty.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://github.com/LithosCarbon/Monty.jl/actions/workflows/documentation.yml)

Monty is a high-performance package for simulating geochemical data produced in experiments, field trials, and commercial deployments in enhanced rock weathering (ERW). See [**the documentation**](https://lithoscarbon.github.io/Monty.jl/) for  explanation, detail, and examples.

## Usage

### Install Julia

To use the `Monty` package, you must have Julia installed. The recommended way to install Julia is with the [juliaup](https://github.com/JuliaLang/juliaup) installer and version manager. You can also follow installation instructions on the [official Julia website](https://julialang.org/downloads/).

### Create an Environment

Before installing the `Monty` package, it's strongly recommended that you create and activate an isolated environment. This can be done by starting Julia and running
```julia
using Pkg
Pkg.activate("my-project")
```
in the directory of your choice, replacing "my-project" with whatever name you want. More information on creating and managing environments can be found [here](https://pkgdocs.julialang.org/v1/environments/#Creating-your-own-environments).

### Install Monty

`Monty` is not in the general Julia registry, so it cannot be installed with `Pkg.add("Monty")`. However, it can easily be installed using the GitHub repository link. With your environment activated, you can run:
```julia
using Pkg
Pkg.add("https://github.com/LithosCarbon/Monty.jl")
```
and the package should automatically be downloaded and installed.

### Load the package

Then, you can use the package by running
```julia
using Monty
```
in a script or at the REPL (again, with your environment activated). This makes all the package's functions and data types available. To see how to get started, look at examples in [the documentation](https://lithoscarbon.github.io/Monty.jl/).
