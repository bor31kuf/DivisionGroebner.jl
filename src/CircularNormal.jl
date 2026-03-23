"""
PolyNome werden mit einem Arrray gespeichert.
Innerhalb der Liste wird ein SIMD Vector bebutzt zur parallelen Operation auf dem Vector.
"""
mutable struct PolyNomCirc3
    Monome::CircularDeque{Vector{Int64}}
    Koeffizienten::CircularDeque{FieldElem}
end



"""
Um das Resultat zu berchnen benutzen wir einen geobucket,für eine schnelle Polynomaddition.
"""
struct geobucketpol3
    Bucket::Vector{PolyNomCirc3}
end


"""
Die Addition bei einem Geobucket
"""
function addgeobucket(B::geobucketpol3,f::PolyNomCirc3)
    #nochmal hinschauen
    i=max(1,ceil(Int,log(4,length(f.Monome))))
    m = length(B.Bucket)
    if i <= m
        B.Bucket[i] =add(B.Bucket[i],f)
        while i <=m && length(B.Bucket[i].Monome) > 4^i
            if i!=m
                B.Bucket[i+1]=add(B.Bucket[i+1],B.Bucket[i])
            else
                push!(B.Bucket,PolyNomCirc3(CircularDeque{Vector{Int64}}(2*4^(m+1)),CircularDeque{FieldElem}(2*4^(m+1))))
                B.Bucket[m+1] = add(B.Bucket[m+1],B.Bucket[m])
            end
            empty!(B.Bucket[i].Monome)
            empty!(B.Bucket[i].Koeffizienten)
            i+=1
        end 
        return B
    end
    for t=m:max(m,i)-1
        push!(B.Bucket, PolyNomCirc3(CircularDeque{Vector{Int64}}(2*4^(t+1)),CircularDeque{FieldElem}(2*4^(t+1))))
    end
    B.Bucket[i] = add(B.Bucket[i],f)
    return B
end


"""
Extrahiert den Leitterm von dem Geobucket
"""
function Leitterm(B::geobucketpol3)
    m= length(B.Bucket)
    j= 0
    while true
        j= 0
        w = true
        for i=1:m
            if isempty(B.Bucket[i].Monome) == false
                if j == 0
                    j=i
                else
                    wt = cmp(first(B.Bucket[i].Monome),first(B.Bucket[j].Monome))
                    if wt==1
                        j=i
                    elseif wt==2
                        if first(B.Bucket[i].Koeffizienten) + first(B.Bucket[j].Koeffizienten)!=0
                            B.Bucket[j].Koeffizienten.buffer[B.Bucket[j].Koeffizienten.first]+=B.Bucket[i].Koeffizienten[1]
                            popfirst!(B.Bucket[i].Koeffizienten)
                            popfirst!(B.Bucket[i].Monome)
                        else
                            popfirst!(B.Bucket[i].Koeffizienten)
                            popfirst!(B.Bucket[i].Monome)
                            popfirst!(B.Bucket[j].Koeffizienten)
                            popfirst!(B.Bucket[j].Monome)
                            w = false
                            break
                        end
                    end
                end
            end
        end
        if j==0 || w== true
            break
        end
    end
    if j== 0
        return PolyNomCirc3(CircularDeque{Vector{Int64}}(1),CircularDeque{FieldElem}(1)) 
    end
    #return
    h = PolyNomCirc3(CircularDeque{Vector{Int64}}(1),CircularDeque{FieldElem}(1))
    push!(h.Monome,popfirst!(B.Bucket[j].Monome))
    push!(h.Koeffizienten,popfirst!(B.Bucket[j].Koeffizienten))
    return h
end

"""
Eine Umwandlung von einem Oscar Polynom in den neuen Polynomtypen.

Unterstützt werden: lex,wdeglex,deglex,degrevlex,wdegrevlex
"""
function PolNeuCirc3O(f;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))+1
    D = PolyNomCirc3(CircularDeque{Vector{Int64}}(L),CircularDeque{FieldElem}(L)) 
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        for i=1:length(A)
            push!(D.Monome,[0,B[i]...])
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdeglex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:deglex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        W = length(gens(parent(f)))*2        
        c = ord.o.matrix
        D = PolyNomCirc3(Vector{Vec{W,Int64}}(),Vector{FieldElem}())
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((c*B[i]...,B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    else
        throw(ArgumentError("Ordnung nicht unterstützt"))
    end
end  


"""
Funktion für den Vergleich von Monomen. 
"""
function cmp(a::Vector{Int64},b::Vector{Int64})
    for i in 1:length(a)
        if a[i] < b[i]
            return 0
        elseif a[i] > b[i]
            return 1
        end
    end
    return 2
end

"""
Funktion zum umwandeln vom neuen Polynomtyp in den Oscar Polynomtypen.
"""
function NeuPolCirc3(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.Monome)
    Builder = MPolyBuildCtx(PolAlg)
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||  typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}
   
        for i=1:k
            push_term!(Builder,f.Koeffizienten[i],reverse(collect(Tuple(f.Monome[i]))[2:end]))
        end
    elseif typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:lex} 
         for i=1:k
            push_term!(Builder,f.Koeffizienten[i],collect(Tuple(f.Monome[i]))[2:end])
        end
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        W = length(gens(PolAlg)) 
        for i=1:k
            push_term!(Builder,f.Koeffizienten[i],collect(Tuple(f.Monome[i]))[W+1:end])
        end
    end
    return finish(Builder)
end




"""
Der eigentliche Divisionsalgortihmus
"""
function DIVCirc3(f::PolyNomCirc3,G::Vector{PolyNomCirc3}) 
    L = length(f.Monome)
    if length(f.Monome)==0
        return f
    end
    f2 = geobucketpol3([PolyNomCirc3(CircularDeque{Vector{Int64}}(8),CircularDeque{FieldElem}(8))])
    f2 =addgeobucket(f2,f)
    LTf2 = Leitterm(f2)
    r = PolyNomCirc3(CircularDeque{Vector{Int64}}(L),CircularDeque{FieldElem}(L))
    D = length(G)
    W = length(G[1].Monome[1])
   
    while length(LTf2.Monome) != 0
        w = false
        for i=1:D
           
            if all(first(LTf2.Monome).>=first(G[i].Monome))
                DIV1 = first(LTf2.Monome)-first(G[i].Monome)
                DIV2 = -first(LTf2.Koeffizienten)/first(G[i].Koeffizienten)

                L2 = length(G[i].Monome)
                A = CircularDeque{Vector{Int64}}(L2-1)
                B = CircularDeque{FieldElem}(L2-1)
                for t=2:L2
                    push!(A,G[i].Monome[t]+DIV1)
                    push!(B,G[i].Koeffizienten[t]*DIV2)
                end
            
                w = true
                
                if length(A)!=0
                    g = PolyNomCirc3(A,B)
                    f2= addgeobucket(f2,g)
                end
                LTf2 = Leitterm(f2)
                break
        
            end
        end
        if w == false
            pushing(r,LTf2)
            LTf2 = Leitterm(f2)
        end
    end
    return r
end

function pushing(r::PolyNomCirc3,LTf2)
    if capacity(r.Koeffizienten) > length(r.Koeffizienten)
        push!(r.Monome,first(LTf2.Monome))
        push!(r.Koeffizienten,first(LTf2.Koeffizienten))
        
        return r
    else
        r21 = CircularDeque{Vec{W,Int64}}(2*capacity(r.Monome))
        r22 = CircularDeque{FieldElem}(2*capacity(r.Monome))
        for i=1:length(r.Koeffizienten)
            push!(r21,r.Monome[i])
            push!(r22,r.Koeffizienten[i]) 
        end
        push!(r.Monome,first(LTf2.Monome))
        push!(r.Koeffizienten,first(LTf2.Koeffizienten))
    end

    return PolyNomCirc(r21,r22)
end


"""
Addition zweier Monome mit Zusatzinfos
"""
function add(f::PolyNomCirc3,g::PolyNomCirc3)
    lf = f.Monome.n
    lg = g.Monome.n
    k= 1
    j= 1
    
    t=0
    tmp = f.Koeffizienten.last +1

    while k <=lf && j <= lg     
    
        x = cmp(f.Monome[k],g.Monome[j])        
        #potentiell aufpassen
        if x == 0
            push!(f.Monome,g.Monome[j])
            push!(f.Koeffizienten,g.Koeffizienten[j])
            j+=1
            
        elseif x==2
            if f.Koeffizienten[k]+g.Koeffizienten[j] != 0
                push!(f.Koeffizienten,f.Koeffizienten[k]+ g.Koeffizienten[j])
                push!(f.Monome,f.Monome[k])
            else
                t-=1
                f.Monome.n +=1
                f.Koeffizienten.n +=1
            end
            k+=1
            j+=1
   
          
        else
            push!(f.Monome,f.Monome[k])
            push!(f.Koeffizienten,f.Koeffizienten[k])
            k+=1
        end
        f.Monome.n -=1
        f.Koeffizienten.n -=1
        t+=1 
    end
    while j <=lg
     
        push!(f.Monome,g.Monome[j])
        push!(f.Koeffizienten,g.Koeffizienten[j])
        j+=1
        t+=1
           
        f.Monome.n -=1
        f.Koeffizienten.n -=1
    end

    while k <=lf
        
        push!(f.Monome,f.Monome[k])
        push!(f.Koeffizienten,f.Koeffizienten[k])
        k+=1
        t+=1
          
        f.Monome.n -=1
        f.Koeffizienten.n -=1
    end 
    f.Koeffizienten.first = ifelse(tmp>f.Koeffizienten.capacity, 1, tmp)
    f.Monome.first = f.Koeffizienten.first
    f.Koeffizienten.n = t
    f.Monome.n = t
  
    return f
end 

"""
Die komplette Divisio
"""
function DIVCircC3(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    f2 = PolNeuCirc3O(f,ord=ord)
    W = length(gens(parent(f)))+1
    G2 = [PolNeuCirc3O(G[i],ord=ord) for i=1:length(G)]
    A = DIVCirc3(f2,G2)
    return NeuPolCirc3(A,parent(f),ord=ord)
end
