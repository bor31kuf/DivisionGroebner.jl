mutable struct QueueElemLex
    Pair::Tuple{Int64,Int64} 
    sugar::Int
    nbits::Int
end


"""
Computes the S-PolyNom
"""

function SPoly(f::PolyNomCircLex{W,T},g::PolyNomCircLex{W,T}) where{W,T}
    kgv =  max(first(f.monoms),first(g.monoms))
    
    mf = kgv-first(f.monoms)
    mg = kgv-first(g.monoms)
    x = addLex(f,g,mf,mg,1/first(f.coefficients),-1/first(g.coefficients)) 
    
    return x
end



"""
The Buchberger Algorithm
"""
function BuchbergerLex(G::Vector{PolyNomCircLex{W,T}},PolAlg) where{W,T}
    L = length(G)
    Queue = BinaryHeap(Base.By(x->(x.sugar,x.nbits)),QueueElemLex[])
    Sugar = [sum(first(G[i].monoms)) for i=1:length(G)]
    for i=2:length(G)
        Queue = QUEUE(G[1:i],Sugar,Queue)
    end
   
    while isempty(Queue) == false
        QueueElem = pop!(Queue)
        Pair = QueueElem.Pair
        sugar = QueueElem.sugar
        
        if Test(G,Pair)==true 
            
            Sij = SPoly(G[Pair[1]],G[Pair[2]])
            if typeof(Sij.monoms) != typeof(G[1].monoms)
                G  = [widen_type(G[i]) for i=1:length(G)] 
            end
            
            S = DIVCircLex(Sij,G) #1 or 2 mmh
            
            println("hi")
            Sk = newPolCircLex(S,PolAlg)
            if length(S.monoms)!=0
                println(length(Queue))
                println(length(G))
                while typeof(S.monoms) != typeof(G[1].monoms)
                    G  = [widen_type(G[i]) for i=1:length(G)] 
              
                end
                t1 = mapreduce(identity,gcd,[numerator(S.coefficients[i]) for i=1:length(S.coefficients)] )
                t2 = mapreduce(identity,gcd,[denominator(S.coefficients[i]) for i=1:length(S.coefficients)] )
                t = QQ(t1,t2)
                S.coefficients.buffer = [divexact!(S.coefficients.buffer[i],t) for i=1:length(S.coefficients.buffer)]
                push!(G,S)  
                push!(Sugar,sugar)
                Queue = QUEUE(G,Sugar,Queue)
               
            end
        end
    end
    return G
end

"""
The gebauer-möller criterias which can be checked before inserting 
"""
function QUEUE(G::Vector{PolyNomCircLex{W,T}},Sugar,Queue) where{W,T}
    
   

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
                tt = maximum( [nbits(G[i].coefficients[j]) for j=1:length(G[i].coefficients)] ) +  maximum( [nbits(G[length(G)].coefficients[j]) for j=1:length(G[length(G)].coefficients)] )
                push!(Queue,QueueElemLex( (i,length(G)) , sugar,tt))
            end
        end
    end
    return Queue


end

"""
Criterais which can be checked after inserting it in the Heap
"""
function Test(G::Vector{PolyNomCircLex{W,T}},Pair) where{W,T}
   
    
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




"""
reduces the groebner basis
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


"""
for a complete reduction 
"""
function interreduceLex(G::Vector{PolyNomCircLex{W,Z}}) where {W,Z}
    L = length(G)
   
    for a=1:L
        w = true
        i = 1
        while i <= length(G)
            b = DIVCircLex(G[i],G[eachindex(G) .!=i])  
                
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





function GroebnerTCircLex(G,PolAlg)
    X= interreduceLex(G)
    X=  BuchbergerLex(G,PolAlg)
    X = reduce_groebner(X)
    X = interreduceLex(X)
    return X
  
end


"""
complete computation with 
"""
function GroebnerCircLex(G;ord=default_ordering(parent(G[1])))
    W =length(gens(parent(G[1])))
    A = [collect(exponents(G[i])) for i=1:length(G)]   
    max_val = maximum(extrema(Iterators.flatten(A))[2])
    Z = minType(max_val)
    T = Vector{PolyNomCircLex{W,Z}}()
    for i=1:length(G)
        push!(T,PolNewCircLex(G[i],Z,ord=ord))
    end
    T=GroebnerTCircLex(T,parent(G[1]))
    
    L =MPolyRingElem[]
    for i=1:length(T)
        push!(L,newPolCircLex(T[i],parent(G[1]),ord=ord))
    end
    L = Oscar.IdealGens(L)
    return L
end




"""
Checks if something is a groebner basis.
"""
function my_isgb(G::Vector{PolyNomCircLex{W,T}}) where{W,T}
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




