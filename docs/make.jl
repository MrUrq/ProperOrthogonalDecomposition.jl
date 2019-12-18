using Documenter, ProperOrthogonalDecomposition

makedocs(;
        modules = [ProperOrthogonalDecomposition],
        format = Documenter.HTML(   assets = ["assets/favicon.ico"],
                                prettyurls = get(ENV, "CI", nothing) == "true"),
        sitename = "ProperOrthogonalDecomposition.jl",
        strict = true,
        clean = true,
        checkdocs = :none,
        pages = Any[
                "Home" => "index.md",
                "Manual" => Any[
                        "man/POD.md",
                        "man/weightedPOD.md",
                        "man/convergence.md",
                ]
        ]
)

deploydocs(
    repo = "github.com/MrUrq/ProperOrthogonalDecomposition.jl.git"
)