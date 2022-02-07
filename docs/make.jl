using Documenter
using SoilMulch

makedocs(
    sitename = "SoilMulch",
    format = Documenter.HTML(),
    modules = [SoilMulch],
    pages=["Home" => "index.md"]

)

deploydocs(;
    repo = "github.com/calabresehydrologylab/SoilMulch.jl",
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
