

"""
Polynoimals are saved in a Circular Deque for fast deleting/inserting at the fron oft the list.
In the Circular Deque there are SIMD-vetors for fast parallel operationson the vector.

The potential weight of the Polynom is also saved.
"""

mutable struct PolyNomCirc{W}
    monoms::CircularDeque{Vec{W,Int64}}
    coefficients::CircularDeque{FieldElem}
end


"""
We are using a geobucket structure, for fast addition of polynomials.

"""
struct geobucketpol1{W}
    bucket::Vector{PolyNomCirc{W}}
end


"""
The addition in a geobucket
"""
function addgeobucket(B::geobucketpol1{W},f::PolyNomCirc,DIV1 = Vec{W,Int64}(ntuple(i-> 0,W)),DIV2 =QQFieldElem(1)) where{W}
    log = cld(64-leading_zeros(length(f.coefficients)),2)
    i=max(1,log)
    m = length(B.bucket)
    if i <= m
        add(B.bucket[i],f,DIV1,DIV2)
        while i <=m && length(B.bucket[i].coefficients) > 4^i
            if i!=m
                add(B.bucket[i+1],B.bucket[i])
                empty!(B.bucket[i].coefficients)
                empty!(B.bucket[i].monoms)
            else
                push!(B.bucket,PolyNomCirc(CircularDeque{Vec{W,Int64}}(2*4^(m+1)),CircularDeque{FieldElem}(2*4^(m+1))))
                add(B.bucket[m+1],B.bucket[m])
                empty!(B.bucket[m].coefficients)
                empty!(B.bucket[m].monoms)
            end
            i+=1
        end
        return B
    end
    for t=m:max(m,i)-1
        push!(B.bucket, PolyNomCirc(CircularDeque{Vec{W,Int64}}(2*4^(t+1)),CircularDeque{FieldElem}(2*4^(t+1))))
    end
    add(B.bucket[i],f,DIV1,DIV2) 
    return B
end


"""
The extraction of the leading term in average log(size of B)
"""
function leading_term(B::geobucketpol1{W}) where{W}
    m= length(B.bucket)
    j= 0
    while true
        j= 0
        w = true
        for i=1:m
            if isempty(B.bucket[i].coefficients) == false
                if j == 0
                    j=i
                else
                    wt = cmp(first(B.bucket[i].monoms),first(B.bucket[j].monoms))
                    if wt==1
                        j=i
                    elseif wt==2
                        add!(B.bucket[j].coefficients.buffer[B.bucket[j].coefficients.first],first(B.bucket[i].coefficients))
                        if iszero(first(B.bucket[j].coefficients))==false
                            popfirst!(B.bucket[i].coefficients)
                            popfirst!(B.bucket[i].monoms)
                        else
                            popfirst!(B.bucket[i].coefficients)
                            popfirst!(B.bucket[i].monoms)
                            popfirst!(B.bucket[j].coefficients)
                            popfirst!(B.bucket[j].monoms)
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
        return 0,0
    end
    return popfirst!(B.bucket[j].monoms),popfirst!(B.bucket[j].coefficients) 
end

"""
A conversion of the Oscar polynomial type to this new one.

supported are: lex,wdeglex,deglex,degrevlex,wdegrevlex
"""
function PolnewCirc(f;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))+1
    
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        D = PolyNomCirc(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}((0,B[i]...)))
            push!(D.coefficients,A[i])
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdeglex}
        D = PolyNomCirc(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
        c = ord.o.weights
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),B[i]...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:deglex}
        D = PolyNomCirc(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),B[i]...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:degrevlex}
        D = PolyNomCirc(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        D = PolyNomCirc(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
        c = ord.o.weights
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.coefficients,A[i])
        end
        return D
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        W = length(gens(parent(f)))*2
        D = PolyNomCirc(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
        c = ord.o.matrix
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}((c*B[i]...,B[i]...)))
            push!(D.coefficients,A[i])
        end
        return D
    else
        throw(ArgumentError("Ordnung nicht unterstützt"))
    end
end  

"""
function for comparing two monomial
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
function for conversion of this new polynomial type to the oscar one
"""
function newPolCirc(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.monoms)
    Builder = MPolyBuildCtx(PolAlg)
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||  typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:k
            push_term!(Builder,f.coefficients[i],reverse(collect(Tuple(f.monoms[i]))[2:end]))
        end
    elseif typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:lex}               
        for i=1:k
            push_term!(Builder,f.coefficients[i],collect(Tuple(f.monoms[i]))[2:end])
        end
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        W = length(gens(PolAlg)) 
        for i=1:k
            push_term!(Builder,f.coefficients[i],collect(Tuple(f.monoms[i]))[W+1:end])
        end
    end
    return finish(Builder)
end


"""
The division algrithm
"""
function DIVCirc(f::PolyNomCirc{W},G::Vector{PolyNomCirc{W}}) where W
    L = length(f.coefficients)
    if L==0
        return f
    end
    f2 = geobucketpol1([PolyNomCirc(CircularDeque{Vec{W,Int64}}(8),CircularDeque{FieldElem}(8))])
    f2 =addgeobucket(f2,f)
    LTf2M= first(f.monoms)
    LTf2K =first(f.coefficients)
    r = PolyNomCirc(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L))
    D = length(G)
    DIV2 = first(G[1].coefficients)
    while true

        w = false
        for i=1:D
            if sum(LTf2M>=first(G[i].monoms))==W
                
                DIV1 =LTf2M-first(G[i].monoms)
                div!(DIV2,-LTf2K,first(G[i].coefficients))
                L2 = length(G[i].coefficients)
                w = true
                if L2!=1
                    f2= addgeobucket(f2,G[i],DIV1,DIV2)
                end
            
                LTf2M,LTf2K = leading_term(f2)
                if iszero(LTf2K)
                    return r
                end
                break
        
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K)
            LTf2M,LTf2K = leading_term(f2)
            if iszero(LTf2K)
                if length(r.coefficients) == 0
                    return r
                end
                if first(r.coefficients) == 1
                    return r
                else 
                    tt = first(r.coefficients)
                    for i =1:length(r.coefficients)
                        div!(r.coefficients.buffer[i],tt)
                    end
                    return r
                end
            end
        end
    end
end

"""
is ignoring an element of G
"""
function DIVCirc(f::PolyNomCirc{W},G::Vector{PolyNomCirc{W}},l) where W
    L = length(f.coefficients)
    if L==0
        return f
    end
    f2 = geobucketpol1([PolyNomCirc(CircularDeque{Vec{W,Int64}}(8),CircularDeque{FieldElem}(8))])
    f2 =addgeobucket(f2,f)
    LTf2M= first(f.monoms)
    LTf2K =first(f.coefficients)
    r = PolyNomCirc(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L))
    D = length(G)
    while true
        w = false
        for i=1:D
            if sum(LTf2M>=first(G[i].monoms))==W && l!=i
                
                DIV1 =LTf2M-first(G[i].monoms)
                div!(LTf2K,-first(G[i].coefficients))
                L2 = length(G[i].coefficients)
                w = true
                if L2!=1
                    f2= addgeobucket(f2,G[i],DIV1,LTf2K)
                end
            
                LTf2M,LTf2K = leading_term(f2)
                if LTf2K == 0
                    return r
                end
                break
        
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K)
            LTf2M,LTf2K = leading_term(f2)
            if LTf2K == 0
                if length(r.coefficients) == 0
                    return r
                end
                if first(r.coefficients) == 1
                    return r
                else 
                    tt = first(r.coefficients)
                    for i =1:length(r.coefficients)
                        r.coefficients.buffer[i] /= tt
                    end
                    return r
                end
            end
        end
    end
end

"""
Because we have a CircularDeque which is fixed in size we have sometimes copy it in a bigger CircularDeque
"""
function pushing(r::PolyNomCirc{W},LTf2M,LTf2K) where{W}
    if capacity(r.coefficients) > length(r.coefficients)
        push!(r.monoms,LTf2M)
        push!(r.coefficients,LTf2K)
        return r
    else
        r21 = CircularDeque{Vec{W,Int64}}(2*capacity(r.monoms))
        r22 = CircularDeque{FieldElem}(2*capacity(r.monoms))
        for i=1:length(r.coefficients)
            push!(r21,r.monoms[i])
            push!(r22,r.coefficients[i]) 
        end
        push!(r21,LTf2M)
        push!(r22,LTf2K)
    end

    return PolyNomCirc(r21,r22)
end


"""
Subtraktion mit Zusatzinfos

"""
function Sub1(f::PolyNomCirc{W},g::PolyNomCirc{W},mf,mg,kf,kg) where{W}
    j=1
    k=1
    lg = length(g.coefficients)
    lf = length(f.coefficients)
    A = CircularDeque{Vec{W,Int64}}(lg+lf)
    C = CircularDeque{FieldElem}(lg+lf)
        

    while k <=lf && j <= lg
      
        x = cmp(f.monoms[k]+mf,g.monoms[j]+mg)
        #potentiell aufpassen
        if x == 0
            push!(A,g.monoms[j]+mg)
            push!(C,g.coefficients[j]*kg)
            j+=1
        elseif x==2
            if g.coefficients[j]*kg+f.coefficients[k]*kf != 0
                push!(C,f.coefficients[k]*kf + g.coefficients[j]*kg)
                push!(A,f.monoms[k]+mf)
            end
            k+=1
            j+=1
        else
            push!(A,f.monoms[k]+mf)
            push!(C,f.coefficients[k]*kf)
            k+= 1
        end    
    end
    while j <=lg
        push!(A,g.monoms[j]+mg)
        push!(C,g.coefficients[j]*kg)
        j+=1
    end
    while k <=lf
        push!(A,f.monoms[k]+mf)
        push!(C,f.coefficients[k]*kf)
        k+=1
    end
    

    f2 = PolyNomCirc(A,C)

    return f2
end


"""
Addition zweier monoms mit Zusatzinfos

so ist das eigentliche Addition f+g*(DIV1,DIV2)
DIV1 ist das Monom, DIV2 der Koeffizient
"""
function add(f::PolyNomCirc{W},g::PolyNomCirc{W},DIV1,DIV2)where{W}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k= 1
    j= 2
    A =f.coefficients.last
    t = 0
    while k <=lf && j <= lg
        t+=1
        m = g.monoms[j]+DIV1
        x = cmp(f.monoms[k],g.monoms[j]+DIV1)
        #potentiell aufpassen
        if x == 0
            #println(" ")
            #println(g.coefficients[j]*DIV2)
            push!(f.monoms,g.monoms[j]+DIV1)
            push!(f.coefficients,g.coefficients[j]*DIV2)
            #println(f.coefficients.buffer[f.coefficients.last])
               
            j+=1
        elseif x==2
            D = f.coefficients[k]
            div!(D,DIV2)
            add!(D,g.coefficients[j])
            mul!(D,DIV2)
            if iszero(D) == false
                push!(f.coefficients,D)
                push!(f.monoms,f.monoms[k])
            else
                f.coefficients.n +=1
                f.monoms.n +=1
                t-=1
            end
            k+=1
            j+=1
        else
            push!(f.monoms,f.monoms[k])
            push!(f.coefficients,f.coefficients[k])
            k+=1
        end
        f.coefficients.n -=1
        f.monoms.n -=1   
    end
    while j <=lg
        push!(f.monoms,g.monoms[j]+DIV1)
        push!(f.coefficients,g.coefficients[j]*DIV2)
        f.coefficients.n -=1
        f.monoms.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.monoms,f.monoms[k])
        push!(f.coefficients,f.coefficients[k])
        f.coefficients.n -=1
        f.monoms.n -=1
        t+=1
        k+=1
    end
    if A+1 <= f.coefficients.capacity
        f.monoms.first = A+1
        f.coefficients.first = A+1 
    else
        f.monoms.first = 1
        f.coefficients.first = 1
    end
    f.monoms.n  = t
    f.coefficients.n = t
    
    return 
    
end

"""
Addition zweier Polynome
"""
function add(f::PolyNomCirc{W},g::PolyNomCirc{W})where{W}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k= 1
    j=1
    A =f.coefficients.last
    t = 0
    while k <=lf && j <= lg
        t+=1
        x = cmp(f.monoms[k],g.monoms[j])
 
        if x == 0
            push!(f.monoms,g.monoms[j])
            push!(f.coefficients,g.coefficients[j])
            j+=1
        elseif x==2
            D  = f.coefficients[k]
            add!(D,g.coefficients[j])   
            if iszero(D) == false
                push!(f.coefficients,D)
                push!(f.monoms,f.monoms[k])
            else
                f.coefficients.n +=1
                f.monoms.n +=1
                t-=1
            end
            k+=1
            j+=1
        else
            push!(f.monoms,f.monoms[k])
            push!(f.coefficients,f.coefficients[k])
            k+=1
        end
        f.coefficients.n -=1
        f.monoms.n -=1   
    end
    while j <=lg
        push!(f.monoms,g.monoms[j])
        push!(f.coefficients,g.coefficients[j])
        f.coefficients.n -=1
        f.monoms.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.monoms,f.monoms[k])
        push!(f.coefficients,f.coefficients[k])
        f.coefficients.n -=1
        f.monoms.n -=1
        t+=1
        k+=1
    end
    if A+1 <= f.coefficients.capacity
        f.monoms.first = A+1
        f.coefficients.first = A+1 
    else
        f.monoms.first = 1
        f.coefficients.first = 1
    end
    f.monoms.n  = t
    f.coefficients.n = t
    return 
    
end

"""
Macht die komplette Division mit Umwandlung davor und danach
"""
function DIVCircC(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    f2 = PolnewCirc(f,ord=ord)
    W = length(gens(parent(f)))+1
    G2 = [PolnewCirc(G[i],ord=ord) for i=1:length(G)]
    A = DIVCirc(f2,G2)
    return newPolCirc(A,parent(f),ord=ord)
end




