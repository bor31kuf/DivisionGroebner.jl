module DivisionGroebner

using Oscar, SIMD,DataStructures

include("CircularSIMDNormal.jl")
include("SIMDArrayNorm.jl")
include("GroebnerSIMDNorm.jl")
include("GroebnerArrayNorm.jl")
include("OhneGeobucket.jl")
include("AlteStruktur.jl")


export DIVCirc, DIVArray, PolNeuCirc, PolNeuArray, NeuPolArray, NeuPolCirc, DIVCircC, DIVArrayC, GroebnerCirc, GroebnerArray, Gewicht, DIVArrayOhne, DivRestFam

end



