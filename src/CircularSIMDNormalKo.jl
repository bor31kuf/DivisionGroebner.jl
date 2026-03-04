"""PolyNome werden mit einer CircularDeque gespeichert, für schnelles löschen am Anfang der Liste.
Innerhalb der Liste wird ein SIMD Vector bebutzt zur parallelen Operation auf dem Vector.

Um nicht immer wieder das potentielle Gewicht des PolyNoms zu berchnen wird es mit gespeichert. 
""" 

mutable struct PolyNomCircKo{W}
    Monome::CircularDeque{Vec{W,Int64}}
    Koeffizienten::CircularDeque{ZZRingElem}
    KoKoeff::ZZRingElem
end


"""
Um das Resultat zu berchnen benutzen wir einen geobucket,für eine schnelle Polynomaddition.
"""
struct geobucketpolKo{W}
    Bucket::Vector{PolyNomCircKo{W}}
end

"""
Die Addition bei einem Geobucket
"""
function addgeobucketKo(B::geobucketpolKo{W},f::PolyNomCircKo,DIV1 = Vec{W,Int64}(ntuple(i-> 0,W)),DIV21=1,DIV22=1) where{W}
    #nochmal hinschauen
    log = cld(64-leading_zeros(length(f.Koeffizienten)),2)
    i=max(1,log)
    m = length(B.Bucket)
    if i <= m
        add(B.Bucket[i],f,DIV1,DIV21,DIV22)
        while i <=m && length(B.Bucket[i].Koeffizienten) > 4^i
            if i!=m
                add(B.Bucket[i+1],B.Bucket[i])
                empty!(B.Bucket[i].Koeffizienten)
                empty!(B.Bucket[i].Monome)
                B.Bucket[i].KoKoeff =1
            else
                push!(B.Bucket,PolyNomCircKo(CircularDeque{Vec{W,Int64}}(2*4^(m+1)),CircularDeque{ZZRingElem}(2*4^(m+1)),ZZ(1)))
                add(B.Bucket[m+1],B.Bucket[m])
                empty!(B.Bucket[m].Koeffizienten)
                empty!(B.Bucket[m].Monome)
            end
            i+=1
        end
        return B
    end
    for t=m:max(m,i)-1
        push!(B.Bucket, PolyNomCircKo(CircularDeque{Vec{W,Int64}}(2*4^(t+1)),CircularDeque{ZZRingElem}(2*4^(t+1)),ZZ(1)))
    end
    add(B.Bucket[i],f,DIV1,DIV21,DIV22) 
    return B
end

function LeittermKo(B::geobucketpolKo{W}) where{W}
    m= length(B.Bucket)
    j= 0
    while true
        j= 0
        w = true
        for i=1:m
            if isempty(B.Bucket[i].Koeffizienten) == false
                if j == 0
                    j=i
                else
                    wt = cmp(first(B.Bucket[i].Monome),first(B.Bucket[j].Monome))
                    if wt==1
                        j=i
                    elseif wt==2
                        if first(B.Bucket[i].Koeffizienten)*B.Bucket[j].KoKoeff + first(B.Bucket[j].Koeffizienten).B.Bucket[i].KoKoeff!=0
                            #aufpassen
                            B.Bucket[j].Koeffizienten.buffer[B.Bucket[j].Koeffizienten.first] +=first(B.Bucket[i].Koeffizienten)
                            println("hier doof")
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
        return 0,0,0
    end
    return popfirst!(B.Bucket[j].Monome),popfirst!(B.Bucket[j].Koeffizienten),B.Bucket[j].KoKoeff
end

"""
Eine Umwandlung von einem Oscar Polynom in den neuen Polynomtypen.

Unterstützt werden: lex,wdeglex,deglex,degrevlex,wdegrevlex
"""
function PolNeuCircKo(f;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))+1
    
    #Aufpassen
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        f =gcd(A)
        D = PolyNomCircKo(CircularDeque{Vec{W,Int64}}(L),CircularDeque{ZZRingElem}(L),ZZ(1)) 
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((0,B[i]...)))
            push!(D.Koeffizienten,Int((A[i]*f.den).num))
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdeglex}
        D = PolyNomCircKo(CircularDeque{Vec{W,Int64}}(L),CircularDeque{ZZRingElem}(L)) 
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:deglex}
        D = PolyNomCircKo(CircularDeque{Vec{W,Int64}}(L),CircularDeque{Int64}(L)) 
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:degrevlex}
        D = PolyNomCircKo(CircularDeque{Vec{W,Int64}}(L),CircularDeque{Int64}(L)) 
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        D = PolyNomCircKo(CircularDeque{Vec{W,Int64}}(L),CircularDeque{Int64}(L)) 
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        W = length(gens(parent(f)))*2
        D = PolyNomCircKo(CircularDeque{Vec{W,Int64}}(L),CircularDeque{Int64}(L)) 
        c = ord.o.matrix
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
function cmp(a::Vec{W,Int64},b::Vec{W,Int64}) where{W}
    @inbounds for i in 1:W
        if a[i] != b[i]
            return a[i] > b[i]
        end
    end
    return 2
end


"""
Funktion zum umwandeln vom neuen Polynomtyp in den Oscar Polynomtypen.
"""
function NeuPolCircKo(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.Monome)
    Builder = MPolyBuildCtx(PolAlg)
    #ganz viel aufpassen
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||  typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:k
            push_term!(Builder,Int(f.Koeffizienten[i]),reverse(collect(Tuple(f.Monome[i]))[2:end]))
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
    g = finish(Builder)
    if g == 0
        return g
    end
    B = collect(coefficients(g,ordering=ord))
    t =gcd(B)
    return g/t
end


"""
Der eigentliche Divisionsalgortihmus
"""
function DIVCircKo(f::PolyNomCircKo{W},G::Vector{PolyNomCircKo{W}}) where W
    L = length(f.Koeffizienten)
    if L==0
        return f
    end
  
    f2 = geobucketpolKo([PolyNomCircKo(CircularDeque{Vec{W,Int64}}(8),CircularDeque{ZZRingElem}(8),ZZ(1))])
    f2 =addgeobucketKo(f2,f)
    LTf2M= first(f.Monome)
    LTf2K1 =first(f.Koeffizienten)
    LTf2K2 = f.KoKoeff
    r = PolyNomCircKo(CircularDeque{Vec{W,Int64}}(L),CircularDeque{ZZRingElem}(L),ZZ(1))
    D = length(G)
    DIV2 = first(G[1].Koeffizienten)
    while true
        w = false
        for i=1:D
            if sum(LTf2M>=first(G[i].Monome))==W
                
                DIV1 =LTf2M-first(G[i].Monome)
                DIV21 = LTf2K1*G[i].KoKoeff
                DIV22 = LTf2K2*first(G[i].Koeffizienten)
                ggt = gcd(DIV21,DIV22)
                DIV21 = divexact(DIV21,ggt)
                DIV22 = divexact(DIV22,ggt)
                L2 = length(G[i].Koeffizienten)
                w = true
                if L2!=1
                    f2= addgeobucketKo(f2,G[i],DIV1,DIV21,DIV22)
                
                end
                
                LTf2M,LTf2K1,LTf2K2 = LeittermKo(f2)
                if LTf2K1 == 0
                  
                    return r
                end
                break
        
            end
        end
        if w == false
            #hier
            r= pushing(r,LTf2M,LTf2K1)
            #Falsch. aber ok
            LTf2M,LTf2K1,LTf2K2 = LeittermKo(f2)
            if LTf2K1 == 0
             
                if length(r.Koeffizienten) == 0
                    return r
                end
                if first(r.Koeffizienten) == 1
                    return r
                else 
                    r.KoKoeff = 1
                    return r
                end
            end
        end
    end
end

  
"""
Weil mit einer CircularDeque gearbeitet wird muss beim Einfügen die Größe potentiell verändert werden.
"""
function pushing(r::PolyNomCircKo{W},LTf2M,LTf2K) where{W}
    if capacity(r.Koeffizienten) > length(r.Koeffizienten)
        push!(r.Monome,LTf2M)
        push!(r.Koeffizienten,LTf2K)
        return r
    else
        r21 = CircularDeque{Vec{W,Int64}}(2*capacity(r.Monome))
        r22 = CircularDeque{ZZRingElem}(2*capacity(r.Monome))
        for i=1:length(r.Koeffizienten)
            push!(r21,r.Monome[i])
            push!(r22,r.Koeffizienten[i]) 
        end
        push!(r21,LTf2M)
        push!(r22,LTf2K)
    end

    return PolyNomCircKo(r21,r22,r.KoKoeff)
end



"""
Addition zweier Monome mit Zusatzinfos

so ist das eigentliche Addition f+g*(DIV1,DIV2)
DIV1 ist das Monom, DIV2 der Koeffizient
"""
function add(f::PolyNomCircKo{W},g::PolyNomCircKo{W},DIV1,DIV21,DIV22)where{W}
    lf = length(f.Koeffizienten)
    lg = length(g.Koeffizienten)
    k= 1
    j= 2
    A =f.Koeffizienten.last
    t = 0
    ggt = gcd(g.KoKoeff*DIV22,f.KoKoeff*DIV21)
    ggt1  = divexact(g.KoKoeff*DIV22,ggt)
    ggt2 = divexact(f.KoKoeff*DIV21,ggt)
    while k <=lf && j <= lg
        t+=1
        x = cmp(f.Monome[k],g.Monome[j]+DIV1)
        #potentiell aufpassen
        if x == 0
            push!(f.Monome,g.Monome[j]+DIV1)
            push!(f.Koeffizienten,g.Koeffizienten[j]*ggt2)
            j+=1
        elseif x==2
            if f.Koeffizienten[k]*ggt1+g.Koeffizienten[j]*ggt2 != 0
                push!(f.Koeffizienten,f.Koeffizienten[k]*ggt1+ g.Koeffizienten[j]*ggt2)
                push!(f.Monome,f.Monome[k])
            else
                f.Koeffizienten.n +=1
                f.Monome.n +=1
                t-=1
            end
            k+=1
            j+=1
        else
            push!(f.Monome,f.Monome[k])
            push!(f.Koeffizienten,f.Koeffizienten[k]*ggt1)
            k+=1
        end
        f.Koeffizienten.n -=1
        f.Monome.n -=1   
    end
    while j <=lg
        push!(f.Monome,g.Monome[j]+DIV1)
        push!(f.Koeffizienten,g.Koeffizienten[j]*ggt2)
        f.Koeffizienten.n -=1
        f.Monome.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.Monome,f.Monome[k])
        push!(f.Koeffizienten,f.Koeffizienten[k]*ggt1)
        f.Koeffizienten.n -=1
        f.Monome.n -=1
        t+=1
        k+=1
    end
    if A+1 <= f.Koeffizienten.capacity
        f.Monome.first = A+1
        f.Koeffizienten.first = A+1 
    else
        f.Monome.first = 1
        f.Koeffizienten.first = 1
    end
    f.Monome.n  = t
    f.Koeffizienten.n = t
    f.KoKoeff *= g.KoKoeff*DIV22
    return 
    
end

"""
Addition zweier Polynome
"""
function add(f::PolyNomCircKo{W},g::PolyNomCircKo{W})where{W}
    lf = length(f.Koeffizienten)
    lg = length(g.Koeffizienten)
    k= 1
    j=1
    A =f.Koeffizienten.last
    t = 0
    ggt = gcd(f.KoKoeff,g.KoKoeff)
    ggt1 = divexact(f.KoKoeff,ggt)
    ggt2 = divexact(g.KoKoeff,ggt)
    while k <=lf && j <= lg
        t+=1
        x = cmp(f.Monome[k],g.Monome[j])
 
        if x == 0
            push!(f.Monome,g.Monome[j])
            push!(f.Koeffizienten,g.Koeffizienten[j]*ggt2)
            j+=1
        elseif x==2
            if f.Koeffizienten[k]*ggt1+g.Koeffizienten[j]*ggt2 != 0
                push!(f.Koeffizienten,f.Koeffizienten[k]*ggt2+ g.Koeffizienten[j]*ggt1)
                push!(f.Monome,f.Monome[k])
            else
                f.Koeffizienten.n +=1
                f.Monome.n +=1
                t-=1
            end
            k+=1
            j+=1
        else
            push!(f.Monome,f.Monome[k])
            push!(f.Koeffizienten,f.Koeffizienten[k]*ggt1)
            k+=1
        end
        f.Koeffizienten.n -=1
        f.Monome.n -=1   
    end
    while j <=lg
        push!(f.Monome,g.Monome[j])
        push!(f.Koeffizienten,g.Koeffizienten[j]*ggt2)
        f.Koeffizienten.n -=1
        f.Monome.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.Monome,f.Monome[k])
        push!(f.Koeffizienten,f.Koeffizienten[k]*ggt1)
        f.Koeffizienten.n -=1
        f.Monome.n -=1
        t+=1
        k+=1
    end
    if A+1 <= f.Koeffizienten.capacity
        f.Monome.first = A+1
        f.Koeffizienten.first = A+1 
    else
        f.Monome.first = 1
        f.Koeffizienten.first = 1
    end
    f.KoKoeff *= g.KoKoeff
    f.Monome.n  = t
    f.Koeffizienten.n = t
    return 
    
end

"""
Macht die komplette Division mit Umwandlung davor und danach
"""
function DIVCircCKo(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    f2 = PolNeuCircKo(f,ord=ord)
    W = length(gens(parent(f)))+1
    G2 = [PolNeuCircKo(G[i],ord=ord) for i=1:length(G)]
    A = DIVCircKo(f2,G2)
    return NeuPolCircKo(A,parent(f),ord=ord)
end
