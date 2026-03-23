"""
PolyNome werden mit einem Arrray gespeichert.
Innerhalb der Liste wird ein SIMD Vector bebutzt zur parallelen Operation auf dem Vector.
"""
mutable struct PolyNomArrayO
    Monome::Vector{Vector{Int64}}
    Koeffizienten::Vector{FieldElem}
end



"""
Um das Resultat zu berchnen benutzen wir einen geobucket,für eine schnelle Polynomaddition.
"""
struct geobucketpolO
    Bucket::Vector{PolyNomArrayO}
end


"""
Die Addition bei einem Geobucket
"""
function addgeobucket(B::geobucketpolO,f::PolyNomArrayO)
    #nochmal hinschauen
    i=max(1,ceil(Int,log(4,length(f.Monome))))
    m = length(B.Bucket)
    if i <= m
        f =add(f,B.Bucket[i])
        while i <=m && length(f.Monome) > 4^i
            if i!=m
                if isempty(B.Bucket[i+1].Koeffizienten) == false
                    f=add(f,B.Bucket[i+1])
                end
            else
                push!(B.Bucket,f)
            end
            empty!(B.Bucket[i].Monome)
            empty!(B.Bucket[i].Koeffizienten)
            i+=1
        end
    end
    for t=m:max(m,i)-1
        push!(B.Bucket, PolyNomArrayO(Vector{Vector{Int64}}(),Vector{FieldElem}()))
    end
    B.Bucket[i] = f
    return B
end


"""
Extrahiert den Leitterm von dem Geobucket
"""
function Leitterm(B::geobucketpolO)
    m= length(B.Bucket)
    j= 0
    while true
        j= 0
        w = true
        for i=1:m
            if first(isempty(B.Bucket[i].Monome)) == false
                if j == 0
                    j=i
                else
                    wt = cmp(first(B.Bucket[i].Monome),first(B.Bucket[j].Monome))
                    if wt==1
                        j=i
                    elseif wt==2
                        if first(B.Bucket[i].Koeffizienten) + first(B.Bucket[j].Koeffizienten)!=0
                            B.Bucket[j].Koeffizienten[1]+=B.Bucket[i].Koeffizienten[1]
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
        return PolyNomArrayO(Vector{Vector{Int64}}(),Vector{FieldElem}()) 
    end
    #return
    h = PolyNomArrayO(Vector{Vector{Int64}}(),Vector{FieldElem}())
    push!(h.Monome,popfirst!(B.Bucket[j].Monome))
    push!(h.Koeffizienten,popfirst!(B.Bucket[j].Koeffizienten))
    return h
end

"""
Eine Umwandlung von einem Oscar Polynom in den neuen Polynomtypen.

Unterstützt werden: lex,wdeglex,deglex,degrevlex,wdegrevlex
"""
function PolNeuArrayO(f;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))+1
    D = PolyNomArrayO(Vector{Vector{Int64}}(),Vector{FieldElem}()) 
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        for i=1:length(A)
            push!(D.Monome,[0,B[i]...])
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdeglex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,[sum(c[j]*B[i][j] for j=1:W-1),B[i]...])
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:deglex}
        for i=1:length(A)
            push!(D.Monome,[sum(B[i][j] for j=1:W-1),B[i]...])
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:length(A)
            push!(D.Monome,[sum(B[i][j] for j=1:W-1),reverse(B[i])...])
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,[sum(c[j]*B[i][j] for j=1:W-1),reverse(B[i])...])
            push!(D.Koeffizienten,A[i])
        end
        return D
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        W = length(gens(parent(f)))*2        
        c = ord.o.matrix
        D = PolyNomArrayO(Vector{Vector{Int64}}(),Vector{FieldElem}())
        for i=1:length(A)
            push!(D.Monome,[c*B[i]...,B[i]...])
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
function NeuPolArrayO(f,PolAlg;ord=default_ordering(PolAlg))
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
function DIVArrayO(f::PolyNomArrayO,G::Vector{PolyNomArrayO}) 
    fk = PolyNomArrayO(copy(f.Monome),copy(f.Koeffizienten))
    L = length(f.Monome)
    if length(f.Monome)==0
        return f
    end
    f2 = geobucketpolO([PolyNomArrayO(Vector{Vector{Int64}}(),Vector{FieldElem}())])

    f2 =addgeobucket(f2,fk)
    LTf2 = Leitterm(f2)
    r = PolyNomArrayO(Vector{Vector{Int64}}(),Vector{FieldElem}())
    D = length(G)
    W = length(G[1].Monome[1])
   
    while length(LTf2.Monome) != 0
        w = false
        for i=1:D
           
            if all(first(LTf2.Monome).>=first(G[i].Monome))
            
                DIV1 = first(LTf2.Monome)-first(G[i].Monome)
                DIV2 = -first(LTf2.Koeffizienten)/first(G[i].Koeffizienten)

                L2 = length(G[i].Monome)
                A = Vector{Vector{Int64}}()
                B = Vector{FieldElem}()
                for t=2:L2
                    push!(A,G[i].Monome[t]+DIV1)
                    push!(B,G[i].Koeffizienten[t]*DIV2)
                end
            
                w = true
                
                if length(A)!=0
                    g = PolyNomArrayO(A,B)
                    f2= addgeobucket(f2,g)
                end
                LTf2 = Leitterm(f2)
                break
                empty!(A)
                empty!(B)
              
                return
            end
        end
        if w == false
            push!(r.Monome,LTf2.Monome[1])
            push!(r.Koeffizienten,LTf2.Koeffizienten[1])
            LTf2 = Leitterm(f2)
        end
    end
    return r
end



"""
Addition zweier Monome mit Zusatzinfos
"""
function add(f::PolyNomArrayO,g::PolyNomArrayO)
    lf = length(f.Monome)
    lg = length(g.Monome)
    k= 1
    j= 1
    #f2 = deepcopy(f)
    A = Vector{Vector{Int64}}()
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
    h = PolyNomArrayO(A,C)
    return h
end

"""
Die komplette Divisio
"""
function DIVArrayCO(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    f2 = PolNeuArrayO(f,ord=ord)
    W = length(gens(parent(f)))+1
    G2 = [PolNeuArrayO(G[i],ord=ord) for i=1:length(G)]
    A = DIVArrayO(f2,G2)
    return NeuPolArrayO(A,parent(f),ord=ord)
end
