module DivisionGroebner

using Oscar, SIMD,DataStructures

include("CircularSIMDNormal.jl")
include("SIMDArrayNorm.jl")
include("GroebnerSIMDNorm.jl")
#include("OscarGeo.jl")
include("OhneGeobucket.jl")
include("AlteStruktur.jl")
include("SIMDCircularNormal.jl")
include("CircularNormal.jl")
include("ArrayNorm.jl")
include("CircularSIMDNormalWeight.jl")


export DIVCirc, DIVArray, PolnewCirc, PolNeuArray, NeuPolArray, newPolCirc, DIVCircC, DIVArrayC, GroebnerCirc, Gewicht, DIVArrayOhneGeoC, DIVArrayOhneGeo, DivRestFam, PolNeuCirc2, NeuPolCirc2, DIVCirc2, DIVCircC2,DivRestFam, add, PolyNomCirc, DIVCirc,DIVArrayCO,DIVCircC3, DIVCircC4, DIVCircCW

end



