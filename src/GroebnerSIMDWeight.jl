mutable struct QueueElemWeight
    Pair::Tuple{Int64,Int64} 
    weight::Int
    nbits::Int
end


"""
Berechnet das S-Polynom
"""
function SPoly(f::PolyNomCircWeight{W,T,Y},g::PolyNomCircWeight{W,T,Y},c) where{W,T,Y}
    kgv =  max(first(f.monoms),first(g.monoms))
    w = sum(kgv*c)
    wf = w-first(f.weight)
    wg = w-first(g.weight)
    mf = kgv-first(f.monoms)
    mg = kgv-first(g.monoms)
    x = SubWeight(f,g,mf,mg,1/first(f.coefficients),-1/first(g.coefficients),wf,wg) 
    return x
end



"""
Der Buchberger Algorithmus
"""
function BuchbergerWeight(G::Vector{PolyNomCircWeight{W,QQFieldElem,T,Y}},PolAlg,c) where{W,T,Y}
    L = length(G)
    Queue = BinaryHeap(Base.By(x->(x.weight,x.nbits)),QueueElemWeight[])
    for i=2:length(G)
        Queue = QUEUE(G[1:i],Queue)
    end
    while isempty(Queue) == false
        QueueElem = pop!(Queue)
        Pair = QueueElem.Pair
        
        if Test(G,Pair)==true 
            
            Sij = SPoly(G[Pair[1]],G[Pair[2]],c)
         
            if typeof(Sij.monoms) != typeof(G[1].monoms)
                G  = [widen_type(G[i]) for i=1:length(G)] 
            end
            if typeof(Sij.weight) != typeof(G[1].weight)
                G  = [widen_type2(G[i]) for i=1:length(G)] 
            end
            S = DIVCircWeight2(Sij,G) #1 oder 2 mmh
          
            if length(S.monoms)!=0
                while typeof(S.monoms) != typeof(G[1].monoms)
                    G  = [widen_type(G[i]) for i=1:length(G)] 
              
                end
                t1 = mapreduce(identity,gcd,[numerator(S.coefficients[i]) for i=1:length(S.coefficients)] )
                t2 = mapreduce(identity,gcd,[denominator(S.coefficients[i]) for i=1:length(S.coefficients)] )
                t = QQ(t1,t2)
                S.coefficients.buffer = [divexact!(S.coefficients.buffer[i],t) for i=1:length(S.coefficients.buffer)]
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
function QUEUE(G::Vector{PolyNomCircWeight{W,QQFieldElem,T,Y}},Queue) where{W,T,Y}
    
   

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
                tt = maximum( [nbits(G[i].coefficients[j]) for j=1:length(G[i].coefficients)] ) +  maximum( [nbits(G[length(G)].coefficients[j]) for j=1:length(G[length(G)].coefficients)] )
                push!(Queue,QueueElemWeight( (i,length(G)) , w,tt))
            end
        end
    end
    return Queue


end


function Test(G::Vector{PolyNomCircWeight{W,QQFieldElem,T,Y}},Pair) where{W,T,Y}
   
    
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


function GroebnerTCircWeight(G,PolAlg,c)
    G= interreduceWeight(G)
    X=  BuchbergerWeight(G,PolAlg,c)
    X = reduce_groebner(X)
    X =interreduceWeight(X)
    return X
  
end

"""
üpberprüft ob etwas ne Gröbnerbasis ist.
"""
function my_isgb(G::Vector{PolyNomCircWeight{W,QQFieldElem,T,Y}}) where{W,T,Y}
    t = length(G)
    for i = 1:t-1
        for j=1:t
            if length(DIVCircWeight(SPoly(G[i],G[j]),G).monoms) != 0
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



function interreduceWeight(G::Vector{PolyNomCircWeight{X,QQFieldElem,Y,Z}}) where {X,Y,Z}
    L = length(G)
   
    for a=1:L
        w = true
        i = 1
        while i <= length(G)
            b = DIVCircWeight(G[i],G[eachindex(G) .!=i])  
                
            while typeof(b.monoms) != typeof(G[1].monoms)
                G  = [widen_type(G[i]) for i=1:length(G)]   
            end
            
            if length(b.monoms)==0
                deleteat!(G,i)
                w = false
                break
            else
                t1 = mapreduce(identity,gcd,[numerator(b.coefficients[i]) for i=1:length(b.coefficients)] )
                t2 = mapreduce(identity,gcd,[denominator(b.coefficients[i]) for i=1:length(b.coefficients)] )
                t = QQ(t1,t2)
                b.coefficients.buffer = [divexact!(b.coefficients.buffer[i],t) for i=1:length(b.coefficients.buffer)]
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
function GroebnerCircWeight(G;ord=default_ordering(parent(G[1])))
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
    
    
    T = Vector{PolyNomCircWeight{W,QQFieldElem,Z,Y}}()
    for i=1:length(G)
        push!(T,PolNewCircWeight(G[i],Z,Y,ord))
    end
    c = GewichtWeight(W,ord,Y)
    T=GroebnerTCircWeight(T,parent(G[1]),c)
    L =MPolyRingElem[]
    for i=1:length(T)
        push!(L,newPolCircWeight(T[i],parent(G[1]),ord=ord))
    end
    L = Oscar.IdealGens(L)
    return L
end

"""
Gibt ds Gewicht wieder.
"""
function GewichtWeight(W,ord,Y)
    c = 0
    if typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}
        c = Vec{W,Y}(ntuple(i->1,W))
    elseif typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex} || typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        c = Vec{W,Y}(tuple(ord.o.weights))
    end
    return c
end



