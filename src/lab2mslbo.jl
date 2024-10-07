module lab2mslbo

using Random
using SDDPlab
using DualSDDP
using SDDP, JuMP

include("aux.jl")
include("bounds.jl")
include("mslbo_builder.jl")

export build_mslbo

end
