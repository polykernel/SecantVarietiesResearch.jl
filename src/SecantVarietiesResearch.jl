module SecantVarietiesResearch

using RelocatableFolders

# Directories
const PKG_DIR = normpath(@__DIR__, "..")
const DATA_DIR = @path joinpath(PKG_DIR, "data")

# Imports
include("imports.jl")

# Exports
include("exports.jl")

# Implementation
include("functions.jl")
include("divisor_db.jl")

end
