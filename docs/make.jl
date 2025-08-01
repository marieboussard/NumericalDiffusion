push!(LOAD_PATH,"../src/")

using Documenter, NumericalDiffusion

makedocs(sitename="NumericalDiffusion",
    pages = [
        "index.md",

        "Finite Volume solving" => ["Getting started" => "fv_solving/getting_started.md",
        "Equations" => "fv_solving/equations.md", "Schemes" => "fv_solving/schemes.md"],

        "Diffusion Quantification" => ["Basic usage" => "numdiff/quantif.md",
        "Projection algorithm" => "numdiff/projection_algo/uzawa.md"]
    ]
)