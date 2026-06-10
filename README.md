# DivisionGroebner.jl

Goal is a division which works with all monomial orders.

Therefore are the following programs


1. AlteStruktur.jl 

This is a primitive function, which just implements the division algorithm with the oscar type and functions. 


>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A1.madi")

>>B =load("B1.madi")

>>DivRestFam(B,A,ord)

This takes 1.9s on average and 408mb of storage

2. OhneGeobucket.jl

This ist still the primitive algorithm, but with our own type.

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A1.madi")

>>B =load("B1.madi")

>>DIVArrayOhneGeoC(B,A,ord)

This takes 150ms on average and 205mb of storage

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A2.madi")

>>B =load("B2.madi")

>>DIVArrayOhneGeoC(B,A,ord)

This takes 1.7s on average and 1.61Gb of storage

3. ArrayNorm.jl 

This uses a geobucket as a Data Structure for fast addition and extracting of the leading term.

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A2.madi")

>>B =load("B2.madi")

>>DIVArrayCO(B,A,ord)

This takes 60 ms on average and 43mb of storage 

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A3.madi")

>>B =load("B3.madi")

>>DIVArrayCO(B,A,ord)

This takes 3.5s on average and 1.91Gb of storage

4. CircularNormal.jl

While extracting the leading term we have to delete the first monomial of the polynomial.
So we use a CircularDeque for fast deleting at the front and fast pushing at the back.

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A3.madi")

>>B =load("B3.madi")

>>DIVCircC3(B,A,ord)

This takes 2.5s on average and 1.03Gb of storage

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A4.madi")

>>B =load("B4.madi")

>>DIVCircC3(B,A,ord)

This takes 9.9s on average and 3.32Gb of storage

>>vars = ["x$i" for i in 1:n]

>>PolAlg, var = polynomial_ring(QQ,vars,internal_ordering=:lex)

>>>>ord = lex(PolAlg)

>>A = load("A5.madi")

>>B =load("B5.madi")

>>DIVCircC3(B,A,ord)

This takes 2.26s on average and 970mb of storage

5. SIMDCircularNormal.jl

Because we have a lot of addition of monomials we use SIMD Vectors to store the exponents.

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A4.madi")

>>B =load("B4.madi")

>>DIVCircC4(B,A,ord)

This takes 8.4s on average and 2.67Gb of storage 

>>vars = ["x$i" for i in 1:n]

>>PolAlg, var = polynomial_ring(QQ,vars,internal_ordering=:lex)

>>>>ord = lex(PolAlg)

>>A = load("A5.madi")

>>B =load("B5.madi")

>>DIVCircC4(B,A,ord)

This takes 3.2s on average and 1.05Gb of storage 

We see that we don't have an time or storage advantage when we have 50 variables

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A6.madi")

>>B =load("B6.madi")

>>DIVCircC4(B,A,ord)

This takes on average 70s and 16.57Gb of storage

6. CircularSIMDNormal.jl 

If we have QQ as the base ring we often have a situation like QQFieldElem = QQFieldElem*QQFieldElem, 
which is expensive.
So we added QQFieldElem(Val(:raw)) in Nemo FlintTypes.jl

function QQFieldElem(::Val{:raw})
    z = new(0,0)
    return z
end

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A6.madi")

>>B =load("B6.madi")

>>DIVCircC(B,A,ord)

This takes on average 55s and 5.51Gb of storage

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A6.madi")

>>B =load("B6.madi")

>>divrem(B,A)

This takes on average 26s and 1.2 Gb of storage, so the in build Oscar function is still faster.

7. CircularSIMDNormalWeight.jl/CircularSIMDNormalLex.jl/CircularSIMDNormalMatrix.jl

For Weight and Matrix, they store the weight now seperata.\\
In Weight every weight ordering and deg ordering is possible, in Lex just Lex and in Matrix every Matrix ordering.

For the weight and monoms we now look first which is the smallest possible type (Int16,Int32...) use it
and change it if there is a overflow. So every ordering is possible.

>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:lex)

>>ord = lex(PolAlg)

>>A = load("A6.madi")

>>B =load("B6.madi")

>>DIVCircCLex(B,A,ord)

This takes on average 35s and 220MiB of storage.


>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:deglex)

>>ord = lex(PolAlg)

>>A = load("A6.madi")

>>B =load("B6.madi")

>>DIVCircCWeight(B,A,ord)

This takes on average 8.9s and 120MiB of storage



>>PolAlg, (x1,x2,x3,x4,x5) = polynomial_ring(QQ,["x1","x2","x3","x4","x5"],internal_ordering=:deglex)

>>A = load("A6.madi")

>>B =load("B6.madi")

>>divrem(B,A)

This takes on average 22s and 1Gb of stroage

8. divison.jl

Puts the Division from 7 together, so division(f,G,ord) works for all ord



