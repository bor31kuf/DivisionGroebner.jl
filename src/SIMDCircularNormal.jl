"""
PolyNome werden mit einem Arrray gespeichert.
Innerhalb der Liste wird ein SIMD Vector bebutzt zur parallelen Operation auf dem Vector.
"""
mutable struct PolyNomCirc4{W}
    Monome::CircularDeque{Vec{W,Int64}}
    coefficients::CircularDeque{FieldElem}
end


"""
Um das Resultat zu berchnen benutzen wir einen geobucket,für eine schnelle Polynomaddition.
"""
struct geobucketpol4{W}
    Bucket::Vector{PolyNomCirc4{W}}
end


"""
Die Addition bei einem Geobucket
"""
function addgeobucket(B::geobucketpol4{W},f::PolyNomCirc4{W}) where {W}
    #nochmal hinschauen
    i=max(1,ceil(Int,log(4,length(f.Monome))))
    m = length(B.Bucket)
    if i <= m
        B.Bucket[i] =add(B.Bucket[i],f)
        while i <=m && length(B.Bucket[i].Monome) > 4^i
            if i!=m
                B.Bucket[i+1]=add(B.Bucket[i+1],B.Bucket[i])
            else
                push!(B.Bucket,PolyNomCirc4(CircularDeque{Vec{W,Int64}}(2*4^(m+1)),CircularDeque{FieldElem}(2*4^(m+1))))
                B.Bucket[m+1] = add(B.Bucket[m+1],B.Bucket[m])
            end
            empty!(B.Bucket[i].Monome)
            empty!(B.Bucket[i].coefficients)
            i+=1
        end 
        return B
    end
    for t=m:max(m,i)-1
        push!(B.Bucket, PolyNomCirc4(CircularDeque{Vec{W,Int64}}(2*4^(t+1)),CircularDeque{FieldElem}(2*4^(t+1))))
    end
    B.Bucket[i] = add(B.Bucket[i],f)
    return B
end


"""
Extrahiert den Leitterm von dem Geobucket
"""
function Leitterm(B::geobucketpol4{W}) where {W}
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
                        if first(B.Bucket[i].coefficients) + first(B.Bucket[j].coefficients)!=0
                            B.Bucket[j].coefficients.buffer[B.Bucket[j].coefficients.first]+=B.Bucket[i].coefficients[1]
                            popfirst!(B.Bucket[i].coefficients)
                            popfirst!(B.Bucket[i].Monome)
                        else
                            popfirst!(B.Bucket[i].coefficients)
                            popfirst!(B.Bucket[i].Monome)
                            popfirst!(B.Bucket[j].coefficients)
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
        return PolyNomCirc4(CircularDeque{Vec{W,Int64}}(1),CircularDeque{FieldElem}(1)) 
    end
    #return
    h = PolyNomCirc4(CircularDeque{Vec{W,Int64}}(1),CircularDeque{FieldElem}(1))
    push!(h.Monome,popfirst!(B.Bucket[j].Monome))
    push!(h.coefficients,popfirst!(B.Bucket[j].coefficients))
    return h
end

"""
Eine Umwandlung von einem Oscar Polynom in den neuen Polynomtypen.

Unterstützt werden: lex,wdeglex,deglex,degrevlex,wdegrevlex
"""
function PolNeuCirc4O(f;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))+1
    D = PolyNomCirc4(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((0,B[i]...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdeglex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),B[i]...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:deglex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),B[i]...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        W = length(gens(parent(f)))*2        
        c = ord.o.matrix
        D = PolyNomCirc4(Vector{Vec{W,Int64}}(),Vector{FieldElem}())
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((c*B[i]...,B[i]...)))
            push!(D.coefficients,A[i])
        end
        return D
    else
        throw(ArgumentError("Ordnung nicht unterstützt"))
    end
end  


"""
Funktion für den Vergleich von Monomen. 
"""
function cmp(a::Vec{W,Int64},b::Vec{W,Int64}) where {W}
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
function NeuPolCirc4(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.Monome)
    Builder = MPolyBuildCtx(PolAlg)
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||  typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}
   
        for i=1:k
            push_term!(Builder,f.coefficients[i],reverse(collect(Tuple(f.Monome[i]))[2:end]))
        end
    elseif typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:lex} 
         for i=1:k
            push_term!(Builder,f.coefficients[i],collect(Tuple(f.Monome[i]))[2:end])
        end
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        W = length(gens(PolAlg)) 
        for i=1:k
            push_term!(Builder,f.coefficients[i],collect(Tuple(f.Monome[i]))[W+1:end])
        end
    end
    return finish(Builder)
end




"""
Der eigentliche Divisionsalgortihmus
"""
function DIVCirc4(f::PolyNomCirc4{W},G::Vector{PolyNomCirc4{W}}) where {W} 
    L = length(f.Monome)
    if length(f.Monome)==0
        return f
    end
    f2 = geobucketpol4([PolyNomCirc4(CircularDeque{Vec{W,Int64}}(8),CircularDeque{FieldElem}(8))])
    f2 =addgeobucket(f2,f)
    f2.Bucket[1].coefficients.buffer = copy(f2.Bucket[1].coefficients.buffer)
    f2.Bucket[1].Monome.buffer = copy(f2.Bucket[1].Monome.buffer)
    LTf2 = Leitterm(f2)
    r = PolyNomCirc4(CircularDeque{Vec{W, Int64}}(L),CircularDeque{FieldElem}(L))
    D = length(G)
   
    while length(LTf2.Monome) != 0
        w = false
        for i=1:D
           
            if sum(first(LTf2.Monome)>=first(G[i].Monome))==W
                DIV1 = first(LTf2.Monome)-first(G[i].Monome)
                DIV2 = -first(LTf2.coefficients)/first(G[i].coefficients)

                L2 = length(G[i].Monome)
                A = CircularDeque{Vec{W,Int64}}(L2-1)
                B = CircularDeque{FieldElem}(L2-1)
                for t=2:L2
                    push!(A,G[i].Monome[t]+DIV1)
                    push!(B,G[i].coefficients[t]*DIV2)
                end
            
                w = true
                
                if length(A)!=0
                    g = PolyNomCirc4(A,B)
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

function DIVCirc4(f::PolyNomCirc4{W},G::Vector{PolyNomCirc4{W}},m) where {W} 
    L = length(f.Monome)
    if length(f.Monome)==0
        return f
    end
    f2 = geobucketpol4([PolyNomCirc4(CircularDeque{Vec{W,Int64}}(8),CircularDeque{FieldElem}(8))])
    f2 =addgeobucket(f2,f)
    f2.Bucket[1].coefficients.buffer = copy(f2.Bucket[1].coefficients.buffer)
    f2.Bucket[1].Monome.buffer = copy(f2.Bucket[1].Monome.buffer)
    LTf2 = Leitterm(f2)
    r = PolyNomCirc4(CircularDeque{Vec{W, Int64}}(L),CircularDeque{FieldElem}(L))
    D = length(G)
   
    while length(LTf2.Monome) != 0
        w = false
        for i=1:D
            if i!= m
                if sum(first(LTf2.Monome)>=first(G[i].Monome))==W
                    DIV1 = first(LTf2.Monome)-first(G[i].Monome)
                    DIV2 = -first(LTf2.coefficients)/first(G[i].coefficients)

                    L2 = length(G[i].Monome)
                    A = CircularDeque{Vec{W,Int64}}(L2-1)
                    B = CircularDeque{FieldElem}(L2-1)
                    for t=2:L2
                        push!(A,G[i].Monome[t]+DIV1)
                        push!(B,G[i].coefficients[t]*DIV2)
                    end
                
                    w = true
                    
                    if length(A)!=0
                        g = PolyNomCirc4(A,B)
                        f2= addgeobucket(f2,g)
                    end
                    LTf2 = Leitterm(f2)
                    break
            
                end
            end
        end
        if w == false
            pushing(r,LTf2)
            LTf2 = Leitterm(f2)
        end
    end
    return r
end

function pushing(r::PolyNomCirc4{W},LTf2) where {W}
    if capacity(r.coefficients) > length(r.coefficients)
        push!(r.Monome,first(LTf2.Monome))
        push!(r.coefficients,first(LTf2.coefficients))
        
        return r
    else
        r21 = CircularDeque{Vec{W,Int64}}(2*capacity(r.Monome))
        r22 = CircularDeque{FieldElem}(2*capacity(r.Monome))
        for i=1:length(r.coefficients)
            push!(r21,r.Monome[i])
            push!(r22,r.coefficients[i]) 
        end
        push!(r.Monome,first(LTf2.Monome))
        push!(r.coefficients,first(LTf2.coefficients))
    end

    return PolyNomCirc(r21,r22)
end

"""
Subtration of 2 Polynomials
"""
function Sub1(f::PolyNomCirc4{W},g::PolyNomCirc4{W},mf,mg,kf,kg) where {W}
    j=1
    k=1
    lg = length(g.coefficients)
    lf = length(f.coefficients)
    A = CircularDeque{Vector{Int64}}(lg+lf)
    C = CircularDeque{FieldElem}(lg+lf)
        

    while k <=lf && j <= lg
      
        x = cmp(f.Monome[k]+mf,g.Monome[j]+mg)
        #potentiell aufpassen
        if x == 0
            push!(A,g.Monome[j]+mg)
            push!(C,g.coefficients[j]*kg)
            j+=1
        elseif x==2
            if g.coefficients[j]*kg+f.coefficients[k]*kf != 0
                push!(C,f.coefficients[k]*kf + g.coefficients[j]*kg)
                push!(A,f.Monome[k]+mf)
            end
            k+=1
            j+=1
        else
            push!(A,f.Monome[k]+mf)
            push!(C,f.coefficients[k]*kf)
            k+= 1
        end    
    end
    while j <=lg
        push!(A,g.Monome[j]+mg)
        push!(C,g.coefficients[j]*kg)
        j+=1
    end
    while k <=lf
        push!(A,f.Monome[k]+mf)
        push!(C,f.coefficients[k]*kf)
        k+=1
    end
    

    f2 = PolyNomCirc4(A,C)

    return f2
end

"""
Addition zweier Polynome mit Zusatzinfos
"""
function add(f::PolyNomCirc4{W},g::PolyNomCirc4{W}) where {W}
    lf = f.Monome.n
    lg = g.Monome.n
    k= 1
    j= 1
    
    t=0
    tmp = f.coefficients.last +1

    while k <=lf && j <= lg     
    
        x = cmp(f.Monome[k],g.Monome[j])        
        #potentiell aufpassen
        if x == 0
            push!(f.Monome,g.Monome[j])
            push!(f.coefficients,g.coefficients[j])
            j+=1
            
        elseif x==2
            if f.coefficients[k]+g.coefficients[j] != 0
                push!(f.coefficients,f.coefficients[k]+ g.coefficients[j])
                push!(f.Monome,f.Monome[k])
            else
                t-=1
                f.Monome.n +=1
                f.coefficients.n +=1
            end
            k+=1
            j+=1
   
          
        else
            push!(f.Monome,f.Monome[k])
            push!(f.coefficients,f.coefficients[k])
            k+=1
        end
        f.Monome.n -=1
        f.coefficients.n -=1
        t+=1 
    end
    while j <=lg
     
        push!(f.Monome,g.Monome[j])
        push!(f.coefficients,g.coefficients[j])
        j+=1
        t+=1
           
        f.Monome.n -=1
        f.coefficients.n -=1
    end

    while k <=lf
        
        push!(f.Monome,f.Monome[k])
        push!(f.coefficients,f.coefficients[k])
        k+=1
        t+=1
          
        f.Monome.n -=1
        f.coefficients.n -=1
    end 
    f.coefficients.first = ifelse(tmp>f.coefficients.capacity, 1, tmp)
    f.Monome.first = f.coefficients.first
    f.coefficients.n = t
    f.Monome.n = t
  
    return f
end 

"""
Die komplette Divisio
"""
function DIVCircC4(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    f2 = PolNeuCirc4O(f,ord=ord)
    W = length(gens(parent(f)))+1
    G2 = [PolNeuCirc4O(G[i],ord=ord) for i=1:length(G)]
    A = DIVCirc4(f2,G2)
    return NeuPolCirc4(A,parent(f),ord=ord)
end
