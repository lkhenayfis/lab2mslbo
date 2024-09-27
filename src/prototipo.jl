module lab2mslbo

using Random
using SDDPlab
using DualSDDP
using SDDP, JuMP

include("helpers.jl")
include("mslbo_builder.jl")

export build_mslbo

end
