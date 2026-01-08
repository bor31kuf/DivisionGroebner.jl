#New Structure for Polynomials, but with 
mutable struct PolyNomArray{W}
    Monome::Vector{Vec{W,Int64}}
    Koeffizienten::Vector{FieldElem}
end

#
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
    end
    if typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdeglex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:deglex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
end  

function Vgl(a::Vec{W,Int64},b::Vec{W,Int64}) where{W}
    for i in 1:W
        if a[i] < b[i]
            return 0
        elseif a[i] > b[i]
            return 1
        end
    end
    return 2
end

function NeuPolArray(f::PolyNomArray,PolAlg) 
    a=zero(PolAlg)
    k = length(f.Monome)
    for i=1:k
        a += monomial(PolAlg,collect(Tuple(f.Monome[i]))[2:end])*f.Koeffizienten[i]
    end
    return a
end




function DIVArray(f::PolyNomArray{W},G::Vector{PolyNomArray{W}}) where W
    fk = PolyNomArray(copy(f.Monome),copy(f.Koeffizienten))
    L = length(f.Monome)
    if length(f.Monome)==0
        return f
    end
    f2 = geobucketpol([PolyNomArray(Vector{Vec{W,Int64}}(),Vector{FieldElem}())])

    f2 =addgeobucket(f2,fk)
    LTf2 = Leitterm(f2)
    r = PolyNomArray(Vector{Vec{W,Int64}}(),Vector{FieldElem}())
    D = length(G)
    while length(LTf2.Monome) != 0
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



function Sub1(f::PolyNomArray{W},g::PolyNomArray{W},mf,mg,kf,kg) where{W}
    j=1
    k=1
    
    lg = length(g.Monome)
    lf = length(f.Monome)
    A = Vector{Vec{W,Int64}}()
    C = Vector{FieldElem}()


    while k <=lf && j <= lg
      
        x = Vgl(f.Monome[k]+mf,g.Monome[j]+mg)
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






struct geobucketpol{W}
    Bucket::Vector{PolyNomArray{W}}
end


function addgeobucket(B::geobucketpol{W},f::PolyNomArray) where{W}
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
        push!(B.Bucket, PolyNomArray(Vector{Vec{W,Int64}}(),Vector{FieldElem}()))
    end
    B.Bucket[i] = f
    return B
end


function add(f::PolyNomArray{W},g::PolyNomArray{W})where{W}
    lf = length(f.Monome)
    lg = length(g.Monome)
    k= 1
    j= 1
    #f2 = deepcopy(f)
    A = Vector{Vec{W,Int64}}()
    C = Vector{FieldElem}()
    while k <=lf && j <= lg
        
        x = Vgl(f.Monome[k],g.Monome[j])

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

function Leitterm(B::geobucketpol{W}) where{W}
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
                    wt = Vgl(first(B.Bucket[i].Monome),first(B.Bucket[j].Monome))
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
        return PolyNomArray(Vector{Vec{W,Int64}}(),Vector{FieldElem}()) 
    end
    #return
    h = PolyNomArray(Vector{Vec{W,Int64}}(),Vector{FieldElem}())
    push!(h.Monome,popfirst!(B.Bucket[j].Monome))
    push!(h.Koeffizienten,popfirst!(B.Bucket[j].Koeffizienten))
    return h
end

function DIVArrayC(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    f2 = PolNeuArray(f,ord=ord)
    W = length(gens(parent(f)))+1
    G2 = [PolNeuArray(G[i],ord=ord) for i=1:length(G)]
    A = DIVArray(f2,G2)
    return NeuPolArray(A,parent(f),ord=ord)
end
