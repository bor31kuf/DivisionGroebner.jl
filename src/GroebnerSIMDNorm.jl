"""
Berechnet das S-Polynom
"""
function SPoly(f::PolyNomCircLex{W,T},g::PolyNomCircLex{W,T}) where{W,T}
    kgv =  max(first(f.monoms),first(g.monoms))
    #kgv = Base.setindex(kgv,sum(kgv*c),1)
    
    mf = kgv-first(f.monoms)
    mg = kgv-first(g.monoms)
    x = Sub1(f,g,mf,mg,1/first(f.coefficients),-1/first(g.coefficients)) 

    return x
end


cmp2(a,b) = cmp(first(a.monoms),first(b.monoms)) == 0 ? true : false

"""
Der Buchberger Algorithmus
"""
function Buchberger3(G::Vector{PolyNomCircLex{W,QQFieldElem,T}},PolAlg) where{W,T}
    L = length(G)
    A = Bool[true for i=1:length(G)]
    Queue = pairs(L)
    Bits = trues(Int(L*(L-1)/2))
    k = 1
    A  = trues(length(G))
    z =length(G)+1
    while k <= length(Bits)
        if Bits[k]
        
            a = Queue[k][1]
            b = Queue[k][2]
            if A[a] && A[b]
                Sij = SPoly(G[a],G[b])
                println("Hi")
                S = DIVCircLex(Sij,G)
                if length(S.monoms)!=0  
                    println(length(G))
                    t = mapreduce(identity,gcd,S.coefficients.buffer[1:S.coefficients.n])
                    S.coefficients.buffer = [divexact!(S.coefficients.buffer[i],t) for i=1:length(S.coefficients.buffer)]
                    push!(G,S) 
                    
                    
                    for i=(z+1):length(G)-1
                         
                        if all(first(S.monoms)<=first(G[i].monoms))    
                            A[i] == false
                        end
                    end
                    push!(A,true)
                    Queue, Bits = QUEUE(G,Queue,Bits,k)
                end
            end
        end
        k+=1
    end
    return G
end

"""
Zeigt auf welche Polynom paare überhaupt in Betracht kommen. 
"""
function QUEUE(G::Vector{PolyNomCircLex{W,QQFieldElem,T}},Pairs,Bits,k) where{W,T}
    
    h = G[length(G)]
    c = length(Bits)
    for i=k+1:length(Bits)
        if Bits[i]
            f=G[Pairs[i][1]]
            g=G[Pairs[i][2]]
            r  = max(first(f.monoms),first(g.monoms))
            w1 = first(h.monoms)<=r
            #w1 = Base.setindex(w1,false,1)
            w2 = max(first(h.monoms),first(f.monoms)) == r
            #w2 = Base.setindex(w2,false,1)
            w3 = max(first(h.monoms),first(g.monoms)) == r
            if sum(w1) == W && sum(w2) !=W && sum(w3) != W #mmmh W-1
                Bits[i] == false
            end
        end
    end

   for i=1:length(G)-1
        push!(Pairs,(i,length(G)))
        w = max(first(G[i].monoms),first(G[length(G)].monoms)) == first(G[i].monoms)+first(G[length(G)].monoms)
        if sum(w) == W
            push!(Bits,false)
        else
            push!(Bits,true)
        end
    end

   for i=1:length(G)-1
        if Bits[c+i]
            for j=i+1:length(G)-1
                if Bits[c+j]
                    r1 = max(first(G[length(G)].monoms),first(G[i].monoms))
                    r2 = max(first(G[length(G)].monoms),first(G[j].monoms))
                    w1 = r1 >= r2
                    w2 = r1 < r2
                    if sum(w1)==W
                        Bits[c+i] =false
                        break
                    elseif sum(w2) == W
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


function GroebnerTCirc(G,PolAlg)
    X=  Buchberger3(G,PolAlg)
    println(my_isgb(G))
    #X = reduce_groebner(X)
    return X
  
end

"""
üpberprüft ob etwas ne Gröbnerbasis ist.
"""
function my_isgb(G::Vector{PolyNomCircLex{W,QQFieldElem,T}}) where{W,T}
    t = length(G)
    for i = 1:t-1
        for j=1:t
            if length(DIVCircLex(SPoly(G[i],G[j]),G).monoms) != 0
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
            a = DIVCircLex(G[i],G[eachindex(G) .!=i])
            if length(a.monoms)==0
                deleteat!(G,i)
                w = false
                break
            else
                w = false
                G[i]=a
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
    W =length(gens(parent(G[1])))
    
    T = Vector{PolyNomCircLex{W,QQFieldElem,Int32}}()
    for i=1:length(G)
        push!(T,PolNewCircLex(G[i],Int32,ord=ord))
    end
    T=GroebnerTCirc(T,parent(G[1]))
    
    L =MPolyRingElem[]
    for i=1:length(T)
        push!(L,newPolCircLex(T[i],parent(G[1]),ord=ord))
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
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        c = ord.o.matrix
    end
    return c
end



