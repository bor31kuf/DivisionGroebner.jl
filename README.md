# DivisionGroebner.jl

Ziel ist eine Division die mit allen Ordnungen funktioniert.

Dabei gibt es in dem Paket folgende
CircularSIMDNormal.jl:l
kann alle Ordnungen benutzen, benutzt geobucket und SIMD Vektoren in einer CircularDeque zum Speichern
PolNeuCirc - wandelt Oscar Polynom in neue Struktur
NeuPolCirc - wandel neue Struktur in Oscar Polynom
DIVCirc - Polynomdivision mit neuer Struktur als Eingabe
DIVCircC -Polynomdivision mit Oscar Polynom als Eingabe,rechnet mit neuer Struktur

SIMDArrayNorm.jl
kann alle Ordnungen benutzen, benutzt geobucket und SIMD Vektoren in einem Array zum Speichern
PolNeuArray - wandelt Oscar Polynom in neue Struktur
NeuPolArray - wandel neue Struktur in Oscar Polynom
DIVArray - Polynomdivision mit neuer Struktur als Eingabe
DIVArrayC -Polynomdivision mit Oscar Polynom als Eingabe,rechnet mit neuer Struktur

OhneGeobucket.jl
kann alle Ordnungen benutzen, benutzt keinen geobucket und SIMD Vektoren in einem Array zum Speichern
PolNeuArray - wandelt Oscar Polynom in neue Struktur
NeuPolArray - wandel neue Struktur in Oscar Polynom
DIVArrayOhneGeo - Polynomdivision mit neuer Struktur als Eingabe
DIVArrayOhneGeoC -Polynomdivision mit Oscar Polynom als Eingabe,rechnet mit neuer Struktur

AlteStruktur.jl
kann alle Ordnungen benutzen, benutzt die Oscar Struktur 
DivRestFam

CircularNormal.jl
kann alle Ordnungen benutzen, benutzt geobucket und normale Vektoren in einer CircularDeque zum Speichern
PolNeuCirc - wandelt Oscar Polynom in neue Struktur
NeuPolCirc - wandel neue Struktur in Oscar Polynom
DIVCirc - Polynomdivision mit neuer Struktur als Eingabe
DIVCircC -Polynomdivision mit Oscar Polynom als Eingabe,rechnet mit neuer Struktur


((CircularSIMDNormal2.jl -))

GroebnerBasen Algorithmus für Gewichtsordnungen (und MatrixOrdnungen)
GroebnerSIMDnorm.jl
GroebnerArraynorm.jl




"""
PolyNomRing, (x1,x2,x3) = polynomial_ring(QQ,["x1","x2","x3"],internal_ordering=:lex)
ord = lex(PolyNomRing)

f = x1^2+x2
G = [x1*x3,x1+x2]

divrem(f,G)
>>(QQMPolyRingElem[0, x1 - x2], x2^2 + x2)

DIVCircC(f,G,ord)
>>x2^2 + x2

DIVArrayC(f,G,ord)
>>x2^2 + x2

DIVArrayOhneGeoC(f,G,ord)
>>x2^2 + x2

DivRestFam(f,G,ord)
>>x2^2 + x2

DIVCircC2(f,G,ord)
>>x2^2 + x2
"""

"""
PolyNomRing, (x1,x2,x3) = polynomial_ring(QQ,["x1","x2","x3"])
ord = matrix_ordering(PolyNomRing,Matrix([1 2 4;2 3 4; 3 4 5]))

f = x1^2+x2
G = [x1*x3,x1+x2]

DIVCircC(f,G,ord)
>>2*x1^2 - 2*x1

DIVArrayC(f,G,ord)
>>2*x1^2 - 2*x1

DIVArrayOhneGeoC(f,G,ord)
>>2*x1^2 - 2*x1

DivRestFam(f,G,ord)
>>2*x1^2 - 2*x1

DIVCircC2(f,G,ord)
>>x1^2 - x1
"""

"""
PolyNomRing, (x1,x2,x3) = polynomial_ring(QQ,["x1","x2","x3"])
ord = deglex(PolyNomRing)

f = x1+x2^2

A = PolNeuCirc(f,ord=ord)
>>DivisionGroebner.PolyNomCirc{4}(CircularDeque{Vec{4, Int64}}([<4 x Int64>[2, 0, 2, 0],<4 x Int64>[1, 1, 0, 0]]), CircularDeque{FieldElem}([1,1]))
NeuPolCirc(A,PolAlg,ord=ord)
>>x1 + x2^2

A = PolNeuCirc2(f,ord=ord)
>>DivisionGroebner.PolyNomCirc2(CircularDeque{Vector{Int64}}([[2, 0, 2, 0],[1, 1, 0, 0]]), CircularDeque{FieldElem}([1,1]))
NeuPolCirc2(A,PolAlg,ord=ord)
>>x1 + x2^2

A = PolNeuArray(f,ord=ord)
>>DivisionGroebner.PolyNomArray{4}(Vec{4, Int64}[<4 x Int64>[2, 0, 2, 0], <4 x Int64>[1, 1, 0, 0]], FieldElem[1, 1])
NeuPolAlrray(A,PolAlg,ord=ord)
>>x1 + x2^2

"""
