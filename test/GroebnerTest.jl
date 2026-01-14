using DivisionGroebner

PolAlg, (x1,x2,x3,x4,x5,x6,x7,x8) = polynomial_ring(QQ,[:x1,:x2,:x3,:x4,:x5,:x6,:x7,:x8],internal_ordering=:lex)

ord = lex(PolAlg)

G= [x1+x2+x3+x4+x5,x1*x2+x2*x3+x3*x4+x4*x5,x1*x2*x3+x2*x3*x4+x5*x3*x4,x1*x2*x3*x4+x5*x2*x3*x4,x1*x2*x3*x4*x5-1]
#G = [x1+x2+x3+x4+x5+x6,x1*x2+x2*x3+x3*x4+x4*x5+x5*x6,x1*x2*x3+x2*x3*x4+x3*x4*x5+x4*x5*x6,x1*x2*x3*x4+x2*x3*x4*x5+x3*x4*x5*x6,x1*x2*x3*x4*x5+x2*x3*x4*x5*x6,x1*x2*x3*x4*x5*x6-1]
#G = [x1+x2+x3+x4+x5+x6+x7,x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x6*x7,x1*x2*x3+x2*x3*x4+x3*x4*x5+x4*x5*x6+x5*x6*x7,x1*x2*x3*x4+x2*x3*x4*x5+x3*x4*x5*x6+x4*x5*x6*x7,x1*x2*x3*x4*x5+x2*x3*x4*x5*x6+x3*x4*x5*x6*x7,x1*x2*x3*x4*x5*x6+x2*x3*x4*x5*x6*x7,x1*x2*x3*x4*x5*x6*x7-1]

println("Hi")
function Test1()
    return collect(groebner_basis(ideal(G),complete_reduction=true,ordering=ord))
end
function Test2()
    return GroebnerCirc(G,ord=lex(PolAlg))
end

