using Documenter
# push!(LOAD_PATH, "../src/")
include("../src/BinaryMatrixProducts.jl")
using .BinaryMatrixProducts

makedocs(
    sitename="BMP Library Documentation",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "BMP" => "man_bmp.md",
        ],
        "Reference" => [
            "BMP" => "ref_bmp.md",
        ]
    ],
    format = Documenter.HTML(
        prettyurls = false
    )
)
