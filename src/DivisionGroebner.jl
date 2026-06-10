module DivisionGroebner

using Oscar, SIMD,DataStructures
include("CircularSIMDNormal.jl")
include("OhneGeobucket.jl")
include("ArrayNorm.jl")
include("CircularNormal.jl")
include("SIMDCircularNormal.jl")
include("CircularSIMDNormal.jl")
include("CircularSIMDNormalLex.jl")
include("CircularSIMDNormalWeigth.jl")
include("CircularSIMDNormalMatrix.jl")
include("division.jl")
include("GroebnerSIMDNorm.jl")


export DIVRestFam, DIVArrayOhneGeoC,DIVArrayC0,DIVCircC3,DIVCircC4, DIVCircCLex,division,GroebnerCirc, DIVCircCWeight, DIVCircCMatrix, DIVCircC
end



