using Documenter, RadMod2D

push!(LOAD_PATH,"../src/")

makedocs(
    modules = [RadMod2D],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://github.com/bueschgens/RadMod2D",
    ),
    sitename = "RadMod2D.jl",
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
)