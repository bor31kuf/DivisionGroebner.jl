module DivisionGroebner

using Oscar, SIMD,DataStructures

include("CircularSIMDNormal.jl")
include("SIMDArrayNorm.jl")
export DIVCirc, DIVArray, PolNeuCirc, PolNeuArray, NeuPolArray, NeuPolCirc, DIVCircC

end



