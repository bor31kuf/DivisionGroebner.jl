module DivisionGroebner

using Oscar, SIMD,DataStructures

include("CircularSIMDNormal.jl")
include("SIMDArrayNorm.jl")
export DivCirc, DIVArray, PolNeuCirc, PolNeuArray, NeuPolArray, NeuPolCirc

end



