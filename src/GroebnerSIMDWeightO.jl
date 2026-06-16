mutable struct QueueElemWeightO
    Pair::Tuple{Int64,Int64} 
    weight::Int
end


"""
Berechnet das S-Polynom
"""
function SPoly(f::PolyNomCircWeightO{W,T,Y},g::PolyNomCircWeightO{W,T,Y},c) where{W,T,Y}
    kgv =  max(first(f.monoms),first(g.monoms))
    w = sum(kgv*c)
    wf = w-first(f.weight)
    wg = w-first(g.weight)
    mf = kgv-first(f.monoms)
    mg = kgv-first(g.monoms)
    x = SubWeightO(f,g,mf,mg,1/first(f.coefficients),-1/first(g.coefficients),wf,wg) 
    return x
end



"""
Der Buchberger Algorithmus
"""
function BuchbergerWeightO(G::Vector{PolyNomCircWeightO{W,FieldElem,T,Y}},PolAlg,c) where{W,T,Y}
    L = length(G)
    Queue = BinaryHeap(Base.By(x->(x.weight)),QueueElemWeightO[])
    for i=2:length(G)
        Queue = QUEUE(G[1:i],Queue)
    end
    while isempty(Queue) == false
        QueueElem = pop!(Queue)
        Pair = QueueElem.Pair
        if Test(G,Pair)==true 
            Sij = SPoly(G[Pair[1]],G[Pair[2]],c)
         
            if typeof(Sij.monoms) != typeof(G[1].monoms)
                G  = [widen_type(G[i]) for i=1:lengthf(G)] 
            end
            if typeof(Sij.weight) != typeof(G[1].weight)
                G  = [widen_type2(G[i]) for i=1:length(G)] 
            end

            S = DIVCircWeightO2(Sij,G) #1 oder 2 mmh
            if length(S.monoms)!=0
                while typeof(S.monoms) != typeof(G[1].monoms)
                    G  = [widen_type(G[i]) for i=1:length(G)] 
              
                end
                push!(G,S)
                Queue = QUEUE(G,Queue)
               
            end
        end
    end
    return G
end

"""
Zeigt auf welc
Hihe Polynom paare überhaupt in Betracht kommen. 
"""
function QUEUE(G::Vector{PolyNomCircWeightO{W,FieldElem,T,Y}},Queue) where{W,T,Y}
    
   

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
                w = max(first(G[i].weight),first(G[length(G)].weight))
                push!(Queue,QueueElemWeightO( (i,length(G)) , w))
            end
        end
    end
    return Queue


end


function Test(G::Vector{PolyNomCircWeightO{W,FieldElem,T,Y}},Pair) where{W,T,Y}
   
    
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


function GroebnerTCircWeightO(G,PolAlg,c)
    G= interreduceWeightO(G)
    X=  BuchbergerWeightO(G,PolAlg,c)
    X = reduce_groebner(X)
    X =interreduceWeightO(X)
    return X
  
end

"""
üpberprüft ob etwas ne Gröbnerbasis ist.
"""
function my_isgb(G::Vector{PolyNomCircWeightO{W,FieldElem,T,Y}}) where{W,T,Y}
    t = length(G)
    for i = 1:t-1
        for j=1:t
            if length(DIVCircWeightO(SPoly(G[i],G[j]),G).monoms) != 0
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



function interreduceWeightO(G::Vector{PolyNomCircWeightO{W,FieldElem,Z,Y}}) where{W,Z,Y}
    L = length(G)
   
    for a=1:L
        w = true
        i = 1
        while i <= length(G)
            b = DIVCircWeightO(G[i],G[eachindex(G) .!=i])  
                
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
    end
    return G
end

"""
Komplette funktion zur Berechnung der Gröbnerbasis
"""
function GroebnerCircWeightO(G;ord=default_ordering(parent(G[1])))
    W =length(gens(parent(G[1])))
    A = [collect(exponents(G[i])) for i=1:length(G)]   
    max_val = maximum(extrema(Iterators.flatten(A))[2])
    Z = minType(max_val)
    
    Y = Int8
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex}|| typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} 
        A= [[sum(A[i][j]*ord.o.weights) for j=1:length(A[i])] for i=1:length(A)]
        B = [sum(B[j]*ord.o.weights) for j=1:length(B)] 
        Y= minType(maximum(extrema(A)[2])) 
    else
        A= [[sum(A[i][j]) for j=1:length(A[i])] for i=1:length(A)]
        Y= minType(maximum(extrema(A)[2])) 
    end
    
    
    T = Vector{PolyNomCircWeightO{W,FieldElem,Z,Y}}()
    for i=1:length(G)
        push!(T,PolNewCircWeightO(G[i],Z,Y,ord))
    end
    c = GewichtWeightO(W,ord,Y)
    T=GroebnerTCircWeightO(T,parent(G[1]),c)
    L =MPolyRingElem[]
    for i=1:length(T)
        push!(L,newPolCircWeightO(T[i],parent(G[1]),ord=ord))
    end
    L = Oscar.IdealGens(L)
    return L
end

"""
Gibt ds Gewicht wieder.
"""
function GewichtWeightO(W,ord,Y)
    c = 0
    if typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}
        c = Vec{W,Y}(ntuple(i->1,W))
    elseif typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex} || typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        c = Vec{W,Y}(tuple(ord.o.weights))
    end
    return c
end



