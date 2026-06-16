mutable struct QueueElemLexO
    Pair::Tuple{Int64,Int64} 
    sugar::Int
end


"""
Berechnet das S-Polynom
"""
function SPoly(f::PolyNomCircLexO{W,T},g::PolyNomCircLexO{W,T}) where{W,T}
    kgv =  max(first(f.monoms),first(g.monoms))
    
    mf = kgv-first(f.monoms)
    mg = kgv-first(g.monoms)
    x = SubLexO(f,g,mf,mg,1/first(f.coefficients),-1/first(g.coefficients)) 
    
    return x
end


"""
Der Buchberger Algorithmus
"""
function BuchbergerLexO(G::Vector{PolyNomCircLexO{W,FieldElem,T}},PolAlg) where{W,T}
    L = length(G)
    Queue = BinaryHeap(Base.By(x->(x.sugar)),QueueElemLexO[])
    Sugar = [sum(first(G[i].monoms)) for i=1:length(G)]
    for i=2:length(G)
        Queue = QUEUEO(G[1:i],Sugar,Queue)
    end
    while isempty(Queue) == false
        QueueElemO = pop!(Queue)
        Pair = QueueElemO.Pair
        sugar = QueueElemO.sugar
        if Test(G,Pair)==true 
            Sij = SPoly(G[Pair[1]],G[Pair[2]])
         
            if typeof(Sij.monoms) != typeof(G[1].monoms)
            
                G  = [widen_type(G[i]) for i=1:length(G)] 
                
            end
 
            S = DIVCircLexO(Sij,G)
            if length(S.monoms)!=0
                while typeof(S.monoms) != typeof(G[1].monoms)
                    G  = [widen_type(G[i]) for i=1:length(G)] 
                end
                push!(G,S)  
                push!(Sugar,sugar)
                Queue = QUEUEO(G,Sugar,Queue)
            
            end
        end
    end
    return G
end

"""
Zeigt auf welc
Hihe Polynom paare überhaupt in Betracht kommen. 
"""
function QUEUEO(G::Vector{PolyNomCircLexO{W,FieldElem,T}},Sugar,Queue) where{W,T}
    
   

   for i=1:length(G)-1
      
        w = max(first(G[i].monoms),first(G[length(G)].monoms)) == first(G[i].monoms)+first(G[length(G)].monoms)
        if sum(w) != W
            w2 = true
            for j=1:length(G)-1
                if i!=j
                    if all( (max( first(G[length(G)].monoms),first(G[j].monoms) )-max( first(G[length(G)].monoms),first(G[i].monoms) )) <=0)
                        if sum( max(first(G[length(G)].monoms),first(G[j].monoms))== max(first(G[length(G)].monoms),first(G[i].monoms))  ) == W
                            if i < j
                                w2 = false
                            end
                        else
                            w2 = false
                        end
                    end
                end
            end
            if w2 == true
                sugar = sum(max(first(G[i].monoms),first(G[length(G)].monoms)) ) +max(Sugar[i]-sum(first(G[i].monoms)),Sugar[length(G)]-sum(first(G[length(G)].monoms)))
                push!(Queue,QueueElemLexO((i,length(G)) , sugar))
            end
        end
    end
    return Queue

DIV
end


function Test(G::Vector{PolyNomCircLexO{W,FieldElem,T}},Pair) where{W,T}
   
    for i=1:length(G)
        r  = max(first(G[Pair[1]].monoms),first(G[Pair[2]].monoms))
        w1 = first(G[i].monoms)<=r
        w2 = max(first(G[i].monoms),first(G[Pair[1]].monoms)) == r
        w3 = max(first(G[i].monoms),first(G[Pair[2]].monoms)) == r

        if sum(w1) == W && sum(w2) !=W && sum(w3) != W 
            return false
        end
    end
    

    return true

end


function GroebnerTCircLexO(G,PolAlg)
    G= interreduceLexO(G)
    X=  BuchbergerLexO(G,PolAlg)
    X = reduce_groebner(X)
    X = interreduceLexO(X)
    return X
  
end

"""
üpberprüft ob etwas    @ REPL[3]:1ne Gröbnerbasis ist.
"""
function my_isgb(G::Vector{PolyNomCircLexO{W,FieldElem,T}}) where{W,T}
    t = length(G)
    for i = 1:t-1
        for j=1:t
            if length(DIVCircLexO(SPoly(G[i],G[j]),G).monoms) != 0
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
    L  = length(G)
    A = trues(L)
    for i=1:L
        for j = 1:L
            if all(first(G[i].monoms)-first(G[j].monoms) >=0) && i!=j
                A[i] = false
                break
            end 
        end
    end
    G= G[A]
    return G
end



function interreduceLexO(G::Vector{PolyNomCircLexO{X,FieldElem,Y}}) where{X,Y}
    L = length(G)
   
 
    i = 1
    while i <= length(G)
        b = DIVCircLexO(G[i],G[eachindex(G) .!=i])  
                
        while typeof(b.monoms) != typeof(G[1].monoms)
            G  = [widen_type(G[i]) for i=1:length(G)]   
        end
            
        if length(b.monoms)==0
            deleteat!(G,i)
            w = false
            break
        else
            w = false
            G[i]=b
            i+=1
        end
    end
    return G
end

"""
Komplette funktion zur Berechnung der Gröbnerbasis
"""
function GroebnerCircLexO(G;ord=default_ordering(parent(G[1])))
    W =length(gens(parent(G[1])))
    A = [collect(exponents(G[i])) for i=1:length(G)]   
    max_val = maximum(extrema(Iterators.flatten(A))[2])
    Z = minType(max_val)
    T = Vector{PolyNomCircLexO{W,FieldElem,Z}}()
    for i=1:length(G)
        push!(T,PolNewCircLexO(G[i],Z,ord=ord))
    end
    T=GroebnerTCircLexO(T,parent(G[1]))
    
    L =MPolyRingElem[]
    for i=1:length(T)
        push!(L,newPolCircLexO(T[i],parent(G[1]),ord=ord))
    end
    L = Oscar.IdealGens(L)
    return L
end





