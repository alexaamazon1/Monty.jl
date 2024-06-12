using Documenter
using Monty

makedocs(;
    sitename="Monty",
    authors="Lithos Carbon",
    modules=[Monty],
    format=Documenter.HTML(assets=["assets/favicon.ico"], collapselevel=3),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Mixing Tutorial" => "mixing.md",
            "Simulation Tutorial" => "tutorial.md",
            "Efficient Simulation" => "efficiency.md",
            "Simple Simulation & Analysis" => "simple.md",
            "Big Simulation & Analysis" => "bigsim.md",
            "Containers" => "containers.md",
        ],
        "Reference" => "reference.md",
    ],
    checkdocs=:none,
    warnonly=true,
)
