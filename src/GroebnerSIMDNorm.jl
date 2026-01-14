"""
Berechnet das S-Polynom
"""
function SPoly(f::PolyNomCirc,g::PolyNomCirc,c::Vec{W,Int64}) where{W}
    kgv =  max(first(f.Monome),first(g.Monome))
    kgv = Base.setindex(kgv,sum(kgv*c),1)
    mf = kgv-first(f.Monome)
    mg = kgv-first(g.Monome)
    x = Sub1(f,g,mf,mg,1/first(f.Koeffizienten),-1/first(g.Koeffizienten)) 

    return x
end

"""
Der Buchberger Algorithmus
"""
function Buchberger2(G::Vector{PolyNomCirc{W}},c::Vec{W,Int64}) where{W}
    L = length(G)
    Queue = pairs(L)
    Bits = trues(Int(L*(L-1)/2))
    k = 1
    while k <= length(Bits)
        if Bits[k]
        
            a = Queue[k][1]
            b = Queue[k][2]
            Sij = SPoly(G[a],G[b],c)
            println(k)
            S = DIVCirc(Sij,G)
            println(length(Bits))
            println(length(G))
            if length(S.Monome)!=0
                push!(G,S)
                Queue, Bits = QUEUE(G,Queue,Bits,k)
            end
    
        end
        k+=1
    end
    return G
end

"""
Zeigt auf welche Polynom paare überhaupt in Betracht kommen. 
"""
function QUEUE(G,Pairs,Bits,k,ord)
    #nach caramba.inria.fr/sem-slides/201409111030
    #EDER,Faugere,Martani,Perry,Roune
    #Seminar of the CARAMEL Team in Nancy, France
    #11.9.2014
    PolAlg
    h = G[length(G)]
    c = length(Bits)
    for i=k+1:length(Bits)
        if Bits[i]
            f=G[Pairs[i][1]]
            g=G[Pairs[i][2]]
            r = lcm(leading_monomial(f,ordering=ord),leading_monomial(g,ordering=ord))
            w1 = divides(r,leading_monomial(h,ordering=ord))
            if w1[1] == true && cmp(leading_monomial(h,ordering=ord),leading_monomial(f,ordering=ord)) != 0 && cmp(leading_monomial(h,ordering=ord),leading_monomial(g,ordering=ord)) != 0
                Bits[i] == 0
            end
        end
    end

    for i=1:length(G)-1
        push!(Pairs,(i,length(G)))
        w = cmp(lcm(leading_monomial(G[i],ordering=ord),leading_monomial(G[length(G)],ordering=ord)),leading_monomial(G[i],ordering=ord)*leading_monomial(G[length(G)],ordering=ord))
        if w == 0
            push!(Bits,false)
        else
            push!(Bits,true)
        end
    end

    for i=1:length(G)-1
        if Bits[c+i]
            for j=i+1:length(G)-1
                if Bits[c+j]
                    r1 = lcm(leading_monomial(G[length(G)],ordering=ord),leading_monomial(G[i],ordering=ord))
                    r2 = lcm(leading_monomial(G[length(G)],ordering=ord),leading_monomial(G[j],ordering=ord))

                    w1 = divides(r1,r2)
                    w2 = divides(r2,r1)

                    if w1[1] ==true
                        Bits[c+i] = false
                        break
                    elseif w2[1] ==true
                        Bits[c+j] = false
                    end
                end
            end
        end
    end
    return Pairs, Bits

end

 

function pairs(n::Int)::Vector{Tuple{Int,Int}}
    t = Vector{Tuple{Int,Int}}()
    for i=1:n
        for j=i+1:n
            push!(t,(i,j))
        end
    end
    return t
end


function Groebner(G,c)
    X=  Buchberger2(G,c)

    X = reduce_groebner(X)
    return X
  
end

"""
üpberprüft ob etwas ne Gröbnerbasis ist.
"""
function my_isgb(G::Vector{PolyNomCirc{W}},c) where{W}
    t = length(G)
    for i = 1:t-1
        for j=1:t
            if length(DIVCirc(SPoly(G[i],G[j],c),G).Monome) != 0
                return false
            end
        end
    end
    return true
end


"""
reduziert die Gröbnerbasis
"""
function reduce_groebner(G)
    i = 1
    w =false
    L  = length(G)
    for a=1:L
        w = true
        while i <= length(G)
            G2 = deepcopy(G)
            deleteat!(G2,i)
            a = DIVCirc(G[i],G2)
            if length(a.Monome)==0
                deleteat!(G,i)
                w = false
                break
            elseif G[i].Koeffizienten != a.Koeffizienten
                w = false
                G[i]=a
                i+=1
            else
                i+=1
            end
        end
    end
    return G
end

"""
Komplette funktion zur Berechnung der Gröbnerbasis
"""
function GroebnerCirc(G;ord=default_ordering(parent(G[1])))
    W =length(gens(parent(G[1])))+1
    
    T = Vector{PolyNomCirc{W}}()
    for i=1:length(G)
        push!(T,PolNeuCirc(G[i],ord=ord))
    end
    c = Gewicht(parent(G[1]),ord)
    T=Groebner(T,c)
    
    L =MPolyRingElem[]
    for i=1:length(T)
        push!(L,NeuPolCirc(T[i],parent(G[1]),ord=ord))
    end
    L = Oscar.IdealGens(L)
    return L
end

"""
Gibt ds Gewicht wieder.
"""
function Gewicht(PolAlg,ord)
    W= length(gens(PolAlg))+1
    c = 0
    if typeof(ord.o) == Oscar.Orderings.SymbOrdering{:lex}
        c = Vec{W,Int64}(ntuple(i->0,W))
    elseif typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex}
        c = Vec{W,Int64}(ntuple(i->1,W))
        Base.setindex(c,0,1)
    elseif typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex} || typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        k = ord.o.weights
        pushfirst!(k,0)
        c = Vec{W,Int64}(k)
    end
    return c
end



