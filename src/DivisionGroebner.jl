module DivisionGroebner

using Oscar, SIMD,DataStructures
include("CircularSIMDNormal.jl")
include("OhneGeobucket.jl")
include("ArrayNorm.jl")
include("CircularNormal.jl")
include("SIMDCircularNormal.jl")
include("CircularSIMDNormal.jl")
include("CircularSIMDNormalLex.jl")
include("CircularSIMDNormalLexO.jl")
include("CircularSIMDNormalWeight.jl")
include("CircularSIMDNormalWeightO.jl")
include("CircularSIMDNormalMatrix.jl")
include("CircularSIMDNormalMatrixO.jl")
include("division.jl")
include("GroebnerSIMDLex.jl")
include("GroebnerSIMDLexO.jl")
include("GroebnerSIMDWeight.jl")
include("GroebnerSIMDWeightO.jl")
include("groebner.jl")

export DIVRestFam, DIVArrayOhneGeoC,DIVArrayC0,DIVCircC3,DIVCircC4, DIVCircCLex, DIVCircCLexO,division,GroebnerCircLex, DIVCircCWeight, DIVCircCMatrix, DIVCircC, PolyNomCircLex,DIVCircCWeightO,DIVCircCMatrixO,DIVCircCLex2, GroebnerCircLexO,GroebnerCircWeight, GroebnerCircWeightO, groebner
end



