"""
PolyNome werden mit einem Arrray gespeichert.
Innerhalb der Liste wird ein SIMD Vector bebutzt zur parallelen Operation auf dem Vector.

Um nicht immer wieder das potentielle Gewicht des PolyNoms zu berchnen wird es mit gespeichert. 
""" 
mutable struct PolyNomArray{W}
    Monome::Vector{Vec{W,Int64}}
    Koeffizienten::Vector{FieldElem}
end


"""
Eine Umwandlung von einem Oscar Polynom in den neuen Polynomtypen.

Unterstützt werden: lex,wdeglex,deglex,degrevlex,wdegrevlex
"""
function PolNeuArray(f;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))+1
    D = PolyNomArray(Vector{Vec{W,Int64}}(),Vector{FieldElem}()) 
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((0,B[i]...)))
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
        D = PolyNomArray(Vector{Vec{W,Int64}}(),Vector{FieldElem}())
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
    for i in 1:W
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
function NeuPolArray(f,PolAlg;ord=default_ordering(PolAlg))
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

lex


"""
Der eigentliche Divisionsalgortihmus
"""
function DIVArrayOhneGeo(f::PolyNomArray{W},G::Vector{PolyNomArray{W}}) where W
    f2 = PolyNomArray(copy(f.Monome),copy(f.Koeffizienten))
    L = length(f.Monome)
    if length(f.Monome)==0
        return f
    end
    LTf2 = PolyNomArray([popfirst!(f2.Monome)],FieldElem[popfirst!(f2.Koeffizienten)])
    r = PolyNomArray(Vector{Vec{W,Int64}}(),Vector{FieldElem}())
    D = length(G)
    while true
        w = false
        for i=1:D
            if sum((first(LTf2.Monome)>=first(G[i].Monome)))==W
               
                DIV1 = first(LTf2.Monome)-first(G[i].Monome)
                DIV2 = -first(LTf2.Koeffizienten)/first(G[i].Koeffizienten)
                

                L2 = length(G[i].Monome)
                A = Vector{Vec{W,Int64}}()
                B = Vector{FieldElem}()
                for t=2:L2
                    push!(A,G[i].Monome[t]+DIV1)
                    push!(B,G[i].Koeffizienten[t]*DIV2)
                end
            
                w = true
                
                if length(A)!=0
                    g = PolyNomArray(A,B)
                    f2 = add(f2,g)
                end
                if length(f2.Monome) == 0
                    return r
                end
                LTf2 =PolyNomArray([popfirst!(f2.Monome)],FieldElem[popfirst!(f2.Koeffizienten)])
        
                break
            end
        end
        if w == false
            push!(r.Monome,LTf2.Monome[1])
            push!(r.Koeffizienten,LTf2.Koeffizienten[1])
            if length(f2.Monome)==0
                return r
            end
            LTf2 = PolyNomArray([popfirst!(f2.Monome)],FieldElem[popfirst!(f2.Koeffizienten)])
            
        end
    end
    return r
end



"""
Subtraktion mit Zusatzinfos

"""
function Sub1(f::PolyNomArray{W},g::PolyNomArray{W},mf,mg,kf,kg) where{W}
    j=1
    k=1
    
    lg = length(g.Monome)
    lf = length(f.Monome)
    A = Vector{Vec{W,Int64}}()
    C = Vector{FieldElem}()


    while k <=lf && j <= lg
      
        x = cmp(f.Monome[k]+mf,g.Monome[j]+mg)
        #potentiell aufpassen
        if x == 0
            push!(A,g.Monome[j]+mg)
            push!(C,g.Koeffizienten[j]*kg)
            j+=1
        elseif x==2
            if g.Koeffizienten[j]*kg+f.Koeffizienten[k]*kf != 0
                push!(C,f.Koeffizienten[k]*kf + g.Koeffizienten[j]*kg)
                push!(A,f.Monome[k]+mf)
            end
            k+=1
            j+=1
        else
            push!(A,f.Monome[k]+mf)
            push!(C,f.Koeffizienten[k]*kf)
            k+= 1
        end    
    end
    while j <=lg
        push!(A,g.Monome[j]+mg)
        push!(C,g.Koeffizienten[j]*kg)
        j+=1
    end
    while k <=lf
        push!(A,f.Monome[k]+mf)
        push!(C,f.Koeffizienten[k]*kf)
        k+=1
    end


    f2 = PolyNomArray(A,C)

    return f2
end

"""
Addition zweier Monome ohne Zusatzinfos
"""
function add(f::PolyNomArray{W},g::PolyNomArray{W})where{W}
    lf = length(f.Monome)
    lg = length(g.Monome)
    k= 1
    j= 1
    A = Vector{Vec{W,Int64}}()
    C = Vector{FieldElem}()
    while k <=lf && j <= lg
        
        x = cmp(f.Monome[k],g.Monome[j])

        #potentiell aufpassen
        if x == 0
            push!(A,g.Monome[j])
            push!(C,g.Koeffizienten[j])
            j+=1
        elseif x==2
            if f.Koeffizienten[k]+g.Koeffizienten[j] != 0
                push!(C,f.Koeffizienten[k]+ g.Koeffizienten[j])
                push!(A,f.Monome[k])
            end
            k+=1
            j+=1
        else
            push!(A,f.Monome[k])
            push!(C,f.Koeffizienten[k])
            k+=1
        end   
    end
    while j <=lg
        push!(A,g.Monome[j])
        push!(C,g.Koeffizienten[j])
        j+=1
    end

    while k <=lf
        push!(A,f.Monome[k])
        push!(C,f.Koeffizienten[k])
        k+=1
    end
    h = PolyNomArray(A,C)
    return h
end

"""
Die komplette Divisio
"""
function DIVArrayOhneGeoC(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    f2 = PolNeuArray(f,ord=ord)
    W = length(gens(parent(f)))+1
    G2 = [PolNeuArray(G[i],ord=ord) for i=1:length(G)]
    A = DIVArrayOhneGeo(f2,G2)
    return NeuPolArray(A,parent(f),ord=ord)
end
