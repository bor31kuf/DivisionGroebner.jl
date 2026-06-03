module DivisionGroebner

using Oscar, SIMD,DataStructures
include("CircularSIMDNormal.jl")
include("SIMDArrayNorm.jl")
include("GroebnerSIMDNorm.jl")
include("OhneGeobucket.jl")
include("AlteStruktur.jl")
include("SIMDCircularNormal.jl")
include("CircularNormal.jl")
include("ArrayNorm.jl")
include("CircularSIMDNormalWeight.jl")
include("CircularSIMDNormalLex.jl")


export GroebnerCirc, DIVCircC12, PolyNomCirc12

end



