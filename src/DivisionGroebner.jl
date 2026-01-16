module DivisionGroebner

using Oscar, SIMD,DataStructures

include("CircularSIMDNormal.jl")
include("SIMDArrayNorm.jl")
include("GroebnerSIMDNorm.jl")
#include("GroebnerArrayNorm.jl")
include("OhneGeobucket.jl")
include("AlteStruktur.jl")
include("CircularNormal.jl")


export DIVCirc, DIVArray, PolNeuCirc, PolNeuArray, NeuPolArray, NeuPolCirc, DIVCircC, DIVArrayC, GroebnerCirc, Gewicht, DIVArrayOhneGeoC, DIVArrayOhneGeo, DivRestFam, PolNeuCirc2, NeuPolCirc2,DIVCirc2,DIVCircC2,DivRestFam

end



