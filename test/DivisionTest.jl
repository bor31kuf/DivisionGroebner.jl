using DivisionGroebner

#Polynomring erstellen

n=5 #Anzahl der Variablen
A1= 20 #Anzahl der Polynome pro Division
A2 = 150 #Anzahl der Divisionen pro Testset
A3 = 12 #Größe des Polynoms durch das dividiert wird durch Exponenten
vars = ["x$i" for i in 1:n]

PolAlg, gens = polynomial_ring(QQ,vars,internal_ordering=:lex)

ord = lex(PolAlg) #Ordnung des Polynomrings
c = Gewicht(PolAlg,ord)


#erstellt die Polynome
A = [[(rand(PolAlg,-1:2,0:4,-10:6)) for j=1:A1] for i=1:A2]
B = []
while length(B) != A2
    m= rand(PolAlg,-1:8,3:11,1:30)^A3
    if m!= []
        push!(B,m)
    end
end




#erstellt die Polynome in neuer 
Ane = [[PolNeuCirc(A[i][j],ord=ord) for j=1:A1] for i = 1:A2]
Bne = [PolNeuCirc(B[i],ord=ord) for i= 1:60]
Ane2 = [[PolNeuArray(A[i][j],ord=ord) for j=1:A1] for i = 1:A2]
Bne2 = [PolNeuCirc(B[i],ord=ord) for i= 1:60]


#Sorgt dafür das wir nicht lauter Polynome haben die Null sind
for i=1:60
    j = 1
    while j<= length(Ane[i])
        while length(Ane[i][j].Monome)==0
            A[i][j] = rand(PolAlg,-1:2,0:4,-10:6)
            Ane[i][j] = PolNeuCirc(A[i][j],ord=ord)
            Ane1[i][j] = PolNeuArray(A[i][j],ord=ord)
        end
        j+=1
    end
end

function TestKorrektheit()
    for i=1:length(Ane)
        if DIVCirC(B[i],A[i])!= divrem(B[i],A[i]) && DIVArrayC(B[i],A[i]) != divrem(B[i],A[i])
            println("Funktioniert nicht")
            return B[i], A[i]
        end
    end
    return "Korrekt"
end

function TestDIVCircC()
    for i=1:length(Ane)
        DIVCircC(B[i],A[i])
    end
    return
end 

function TestDIVCirc()
    for i=1:length(Ane)
        DIVCirc(Bne[i],Ane[i])
    end
    return 
end 


function TestNormal()
    for i=1:length(Ane)
        divrem(B[i],A[i])
    end
    return 
end

function TestArrayC()
    for i=1:length(Ane2)
        DIVArrayC(B[i],A[i])
    end
    return
end

function TestArray()
    for i=1:length(Ane2)
        DIVArray(Bne2[i],Ane2[i])
    end
    return
end


