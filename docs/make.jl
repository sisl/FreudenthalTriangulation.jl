push!(LOAD_PATH, "../src/")

using Documenter, FreudenthalTriangulation

makedocs(
    modules = [FreudenthalTriangulation],
    format = Documenter.HTML(),
    sitename = "FreudenthalTriangulation.jl",
    pages = [
        ##############################################
        ## MAKE SURE TO SYNC WITH docs/src/index.md ##
        ##############################################
        "Table of Contents" => [
            "index.md",
            "install.md",
            "usage.md",
            "concepts.md"
           ]
    ]
)

deploydocs(
    repo = "github.com/SidhartK/TestRepo2.jl.git",
)
