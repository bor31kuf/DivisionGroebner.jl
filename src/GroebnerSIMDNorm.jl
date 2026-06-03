"""
Berechnet das S-Polynom
"""
function SPoly(f::PolyNomCirc12{W,T},g::PolyNomCirc12{W,T}) where{W,T}
    kgv =  max(first(f.monoms),first(g.monoms))
    #kgv = Base.setindex(kgv,sum(kgv*c),1)
    
    mf = kgv-first(f.monoms)
    mg = kgv-first(g.monoms)
    x = Sub1(f,g,mf,mg,1/first(f.coefficients),-1/first(g.coefficients)) 

    return x
end

"""
Der Buchberger Algorithmus
"""
function Buchberger3(G::Vector{PolyNomCirc12{W,T}},PolAlg) where{W,T}
    L = length(G)
    Queue = pairs(L)
    Bits = trues(Int(L*(L-1)/2))
    k = 1
    while k <= length(Bits)
        if Bits[k]
        
            a = Queue[k][1]
            b = Queue[k][2]
            Sij = SPoly(G[a],G[b])
      
            S = DIVCirc12(Sij,G)
            
            if length(S.monoms)!=0
                ##
                t1 =gcd([S.coefficients.buffer[i].num for i=1:length(S.coefficients.buffer)])
                t2 =gcd([S.coefficients.buffer[i].den for i=1:length(S.coefficients.buffer)])
                t = QQ(t2,t1)
                S.coefficients.buffer = [mul!(S.coefficients.buffer[i],t) for i=1:length(S.coefficients.buffer)]
                println(S)
                if lenght(G) == 20
                    return
                end
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
function QUEUE(G::Vector{PolyNomCirc12{W,T}},Pairs,Bits,k) where{W,T}
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
    #println(my_isgb(G))
    #X = reduce_groebner(X)
    return X
  
end

"""
üpberprüft ob etwas ne Gröbnerbasis ist.
"""
function my_isgb(G::Vector{PolyNomCirc12{W,T}}) where{W,T}
    t = length(G)
    for i = 1:t-1
        for j=1:t
            if length(DIVCirc12(SPoly(G[i],G[j]),G).monoms) != 0
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
        while i <= length(G)DIVCirc
            a = DIVCirc12(G[i],G,i)
            if length(a.monoms)==0
                deleteat!(G,i)
                w = false
                break
            elseif G[i].coefficients != a.coefficients
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
    W =length(gens(parent(G[1])))
    
    T = Vector{PolyNomCirc12{W,Int32}}()
    for i=1:length(G)
        push!(T,PolnewCirc12(G[i],Int32,ord=ord))
    end
    T=GroebnerTCirc(T,parent(G[1]))
    
    L =MPolyRingElem[]
    for i=1:length(T)
        push!(L,newPolCirc12(T[i],parent(G[1]),ord=ord))
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



