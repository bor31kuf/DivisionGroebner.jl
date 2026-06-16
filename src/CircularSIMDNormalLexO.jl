
"""
Polynoimals are saved in a Circular Deque for fast deleting/inserting at the fron oft the list.
In the Circular Deque there are SIMD-vetors for fast parallel operationson the vector.

The potential weight of the Polynom is also saved.
f"""

mutable struct PolyNomCircLexO{W,T,Z}
    monoms::CircularDeque{Vec{W,Z}}
    coefficients::CircularDeque{T}
end


"""
We are using a geobucket structure, for fast addition of polynomials.

"""
struct geobucketpolLexO{W,Z}
    bucket::Vector{PolyNomCircLexO{W,FieldElem,Z}}
end


"""1
The addition in a geobucket
"""
function addgeobucketLexO(B::geobucketpolLexO{W,Z},f::PolyNomCircLexO{W,FieldElem,Z},DIV1 = Vec{W,Z}(ntuple(i-> 0,W)),DIV2 =one(parent(f.coefficients[1]))) where{W,Z}
    log = cld(64-leading_zeros(length(f.coefficients)),2)
    i=max(1,log)
    m = length(B.bucket)
    if i <= m
        if addLexO(B.bucket[i],f,DIV1,DIV2) == false
            return B,false
        end
        while i <=m && length(B.bucket[i].coefficients) > 4^i
            if i!=m
                addLexO(B.bucket[i+1],B.bucket[i])
                empty!(B.bucket[i].coefficients)
                empty!(B.bucket[i].monoms)
            else
                v = [one(parent(f.coefficients[1])) for z = 1:2*4^(m+1)]
                push!(B.bucket,PolyNomCircLexO(CircularDeque{Vec{W,Z}}(2*4^(m+1)),CircularDeque{FieldElem}(2*4^(m+1))))
                B.bucket[m+1].coefficients.buffer = v
                addLexO(B.bucket[m+1],B.bucket[m])
                empty!(B.bucket[m].coefficients)
                empty!(B.bucket[m].monoms)
            end
            i+=1
        end
        return B,true
    end
    for t=m:max(m,i)-1
        v = [one(parent(f.coefficients[1])) for z = 1:2*4^(t+1)]
        push!(B.bucket, PolyNomCircLexO(CircularDeque{Vec{W,Z}}(2*4^(t+1)),CircularDeque{FieldElem}(2*4^(t+1))))
        B.bucket[t+1].coefficients.buffer = v
    end
    if addLexO(B.bucket[i],f,DIV1,DIV2) == false
        return B,false
    end 
    return B,true
end


"""
The extraction of the leading term in average log(size of B)
"""
function leading_term(B::geobucketpolLexO{W,Z},LTf2K) where{W,Z}
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
                            popfirst2!(B.bucket[i].coefficients)
                            popfirst!(B.bucket[i].monoms)
                        else
                            popfirst2!(B.bucket[i].coefficients)
                            popfirst!(B.bucket[i].monoms)
                            popfirst2!(B.bucket[j].coefficients)
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
        return  Vec{W,Z}(ntuple(i-> 0,W)),QQ(0)
    end
    
    Nemo.set!(LTf2K,popfirst3!(B.bucket[j].coefficients))
    return popfirst!(B.bucket[j].monoms),LTf2K
end

function popfirst3!(D)
    v = first(D)
    D.n -= 1
    tmp = D.first + 1
    D.first = ifelse(tmp > D.capacity, 1, tmp)
    v
end

function popfirst2!(D)
    D.n -= 1
    tmp = D.first + 1
    D.first = ifelse(tmp > D.capacity, 1, tmp)
end



"""
A conversion of the Oscar polynomial type to this new one.

supported are: lex,wdeglex,deglex,degrevlex,wdegrevlex
"""
function PolNewCircLexO(f,Z::Type;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))
    D = PolyNomCircLexO(CircularDeque{Vec{W,Z}}(L),CircularDeque{FieldElem}(L)) 
    for i=1:length(A)
        push!(D.monoms,convert(Vec{W, Z}, vload(Vec{W, Int64}, B[i], 1)))
        push!(D.coefficients,A[i])
    end
    return D
end  

"""
function for comparing two monomial
"""
function cmp(a::Vec{W,Z},b::Vec{W,Z}) where{W,Z}
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
function newPolCircLexO(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.monoms)
    Builder = MPolyBuildCtx(PolAlg)
    for i=1:k
            push_term!(Builder,f.coefficients[i],Int64.(collect(Tuple(f.monoms[i]))[1:end]))
        end
    return finish(Builder)
end


function DIVCircLexO(f::PolyNomCircLexO{W,FieldElem,Z},G::Vector{PolyNomCircLexO{W,FieldElem,Z}}) where {W,Z}
    L = length(f.coefficients)
    if L==0
        return f
    end
    f2 = geobucketpolLexO([PolyNomCircLexO(CircularDeque{Vec{W,Z}}(8),CircularDeque{FieldElem}(8))])
    f2.bucket[1].coefficients.buffer = [one(parent(f.coefficients[1])) for z=1:8]
    
    f2 =addgeobucketLexO(f2,f)[1]
    LTf2M= first(f.monoms)
    LTf2K =first(f.coefficients)
    r = PolyNomCircLexO(CircularDeque{Vec{W,Z}}(L),CircularDeque{FieldElem}(L))
    r.coefficients.buffer = [one(parent(f.coefficients[1])) for z=1:L]

    r,w,LTf2M,LTf2K = CircCirc(LTf2M,LTf2K,f2,G,r)
   
    if w == true
        return r
    else
        return widenProblem(LTf2M,LTf2K,f2,G,r)
    end
end

"""
The division algrithm
"""
function CircCirc(LTf2M::Vec{W,Z},LTf2K::FieldElem,f2::geobucketpolLexO{W,Z},G::Vector{PolyNomCircLexO{W,FieldElem,Z}},r::PolyNomCircLexO{W,FieldElem,Z}) where {W,Z}
    D = length(G)
    DIV2 =one(parent(LTf2K))
    while true
       
        w = false
        for i=1:D
            if all(LTf2M>=first(G[i].monoms))
                
                DIV1 =LTf2M-first(G[i].monoms)  
                DIV2 = LTf2K/first(G[i].coefficients)
           
                neg!(DIV2)
             
                L2 = length(G[i].coefficients)
                w = true
                
                if L2!=1
                    
                    f2, w2 =addgeobucketLexO(f2,G[i],DIV1,DIV2)
                    if w2==false
                        return r, false, LTf2M,LTf2K
                    end
                end
            
                LTf2M,LTf2K = leading_term(f2,LTf2K)
                if iszero(LTf2K)
                    return r, true, LTf2M,LTf2K
                end
                break
        
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K)
            LTf2M,LTf2K = leading_term(f2,LTf2K)
            if iszero(LTf2K)
                return r, true, LTf2M,LTf2K
            end
        end
    end
end

"""
is ignoring an element of G
"""




"""
Because we have a CircularDeque which is fixed in size we have sometimes copy it in a bigger CircularDeque
"""
function pushing(r::PolyNomCircLexO{W,FieldElem,Z},LTf2M,LTf2K) where{W,Z}
    if capacity(r.coefficients) > length(r.coefficients)
        push!(r.monoms,LTf2M)        
        Nemo.set!(r.coefficients.buffer[r.monoms.n],LTf2K)
        r.coefficients.last +=1
        r.coefficients.n +=1
   
        return r
    else
        r21 = CircularDeque{Vec{W,Z}}(2*capacity(r.monoms))
        r22 = CircularDeque{FieldElem}(2*capacity(r.monoms))
        r22.buffer = [one(parent(LTf2K)) for z=1:2*capacity(r.monoms)]
        for i=1:length(r.coefficients)
            push!(r21,r.monoms[i])
            push!(r22,r.coefficients[i]) 
        end
        push!(r21,LTf2M)
        
        Nemo.set!(r22.buffer[r21.n],LTf2K)
        r22.last +=1
        r22.n +=1
    end
    
    return PolyNomCircLexO(r21,r22)
end


"""
Addition zweier monoms mit Zusatzinfos

so ist das eigentliche Addition f+g*(DIV1,DIV2)
DIV1 ist das Monom, DIV2 der Koeffizient
"""
function addLexO(f::PolyNomCircLexO{W,FieldElem,Z},g::PolyNomCircLexO{W,FieldElem,Z},DIV1,DIV2)where{W,Z}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k= 1
    j= 2
    A =f.coefficients.last
    t = 0
    t2 = f.coefficients.buffer
    B = f.monoms.first


    for i=2:g.monoms.n
        if any(g.monoms[i]+DIV1  < 0)
            return false
        end
    end
    while k <=lf && j <= lg
        t+=1
        x = cmp(f.monoms[k],g.monoms[j]+DIV1)
        #potentiell aufpassen
        if x == 0
           
            push!(f.monoms,g.monoms[j]+DIV1)
            mul!(t2[f.monoms.last],g.coefficients[j],DIV2)
               
            j+=1
        elseif x==2
            add!(f.coefficients[k], g.coefficients[j]*DIV2)
            
            if iszero(f.coefficients[k]) == false
                push!(f.monoms,f.monoms[k])
                o = f.coefficients.buffer[f.monoms.last]
                f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
                x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
                f.coefficients.buffer[x] =  o
            else
                f.monoms.n +=1
                t-=1
            end
 
            k+=1
            j+=1
        else
            push!(f.monoms,f.monoms[k])
            o = f.coefficients.buffer[f.monoms.last]
            f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
            x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
            f.coefficients.buffer[x] =  o
            k+=1
        end
        
        f.monoms.n -=1   
    end
    while j <=lg
        push!(f.monoms,g.monoms[j]+DIV1)
        mul!(t2[f.monoms.last],g.coefficients[j],DIV2)
        f.monoms.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.monoms,f.monoms[k])
        o = f.coefficients.buffer[f.monoms.last]
        f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
        x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
        f.coefficients.buffer[x] =  o
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
    f.coefficients.last = f.monoms.last
    f.coefficients.buffer = t2
    return true
    
end


function SubLexO(f::PolyNomCircLexO{W,FieldElem,T},g::PolyNomCircLexO{W,FieldElem,T},DIV1::Vec{W,T},DIV3::Vec{W,T},DIV2::FieldElem,DIV4::FieldElem)where{W,T}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k = 1
    j = 1
   
    T2 = T
    for i=1:g.monoms.n
        if any(g.monoms[i]+DIV3  < 0)
            g= widen_type(g)
            f = widen_type(f)
            T2 = widen(T)
        end
    end
       
    for i=1:f.monoms.n
        if any(f.monoms[i]+DIV1  < 0)
            g= widen_type(g)
            f = widen_type(f)
            T2= widen(T)
        end
    end

    A = CircularDeque{Vec{W,T2}}(lf+lg)
    B = CircularDeque{FieldElem}(lf+lg)
    v =  [one(parent(f.coefficients[1])) for z = 1:(lf+lg)]
    B.buffer = v 
     
    while k <=lf && j <= lg
        x = cmp(f.monoms[k]+DIV1,g.monoms[j]+DIV3)
        #potentiell aufpassen
        if x == 0
            push!(A,g.monoms[j]+DIV3)
            mul!(B.buffer[A.last],g.coefficients[j],DIV4)
               
            j+=1
        elseif x==2
            if iszero(f.coefficients[k]*DIV2+DIV4*g.coefficients[j]) == false
                push!(A,g.monoms[j]+DIV3)
                B.buffer[A.last] = f.coefficients[k]*DIV2+DIV4*g.coefficients[j]     
            end
 
            k+=1
            j+=1
        else
            push!(A,f.monoms[k]+DIV1)
            mul!(B.buffer[A.last],f.coefficients[k],DIV2)
            k+=1
        end
          
    end
    while j <=lg
        push!(A,g.monoms[j]+DIV3)
        mul!(B.buffer[A.last],g.coefficients[j],DIV4)

        j+=1
    end

    while k <=lf
        push!(A,f.monoms[k]+DIV1)
        mul!(B.buffer[A.last],f.coefficients[k],DIV2)
        k+=1
    end
    B.n = A.n
    B.last =A.last
    B.first = A.first
    return PolyNomCircLexO(A,B)
    
end
     

"""
Addition zweier Polynome
"""
function addLexO(f::PolyNomCircLexO{W,FieldElem,Z},g::PolyNomCircLexO{W,FieldElem,Z})where{W,Z}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k= 1
    j=1
    A =f.coefficients.last
    B = f.monoms.first
    t = 0
    while k <=lf && j <= lg
        t+=1
        x = cmp(f.monoms[k],g.monoms[j])
 
        if x == 0
            push!(f.monoms,g.monoms[j])
            o = f.coefficients.buffer[f.monoms.last]
            f.coefficients.buffer[f.monoms.last] = g.coefficients[j]
            x = ifelse(j-1+g.monoms.first> g.monoms.capacity,j-1+g.monoms.first-g.monoms.capacity,j-1+g.monoms.first)
            g.coefficients.buffer[x] =  o
            j+=1
        elseif x==2
            add!(f.coefficients[k],g.coefficients[j])
           
            if iszero(f.coefficients[k]) == false
                push!(f.monoms,f.monoms[k])
                o = f.coefficients.buffer[f.monoms.last]
                f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
                x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
                f.coefficients.buffer[x] =  o
            else
                f.monoms.n +=1
                t-=1
            end
            k+=1
            j+=1
        else
            push!(f.monoms,f.monoms[k])
            o = f.coefficients.buffer[f.monoms.last]
            f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
            x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
            f.coefficients.buffer[x] =  o
            k+=1
        end
        
        f.monoms.n -=1   
    end
    while j <=lg
        push!(f.monoms,g.monoms[j])
        o = f.coefficients.buffer[f.monoms.last]
        f.coefficients.buffer[f.monoms.last] = g.coefficients[j]
        x = ifelse(j-1+g.monoms.first> g.monoms.capacity,j-1+g.monoms.first-g.monoms.capacity,j-1+g.monoms.first)
        g.coefficients.buffer[x] =  o
        f.monoms.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.monoms,f.monoms[k])
        f.monoms.n -=1
        o = f.coefficients.buffer[f.monoms.last]
        f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
        x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
        f.coefficients.buffer[x] =  o
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
    f.coefficients.last = f.monoms.last
    f.coefficients.n = t
    return f
    
end


function widenProblem(LTf2M,LTf2K,f2,G::Vector{PolyNomCircLexO{W,FieldElem,Z}},r) where {W,Z<:Integer}
    LTf2M = convert(Vec{W, widen(Z)}, LTf2M)
    NewVecType = Vec{W, widen(Z)}
    f2 = geobucketpolLexO([widen_type(f2.bucket[i]) for i=1:length(f2.bucket)])
    G  = [widen_type(G[i]) for i=1:length(G)]
    r = widen_type(r)
    r, w, LTf2M, LTf2K = CircCirc2(LTf2M,LTf2K,f2,G,r)
 
    if w == true
        return r
    else
        return widenProblem(LTf2M,LTf2K,f2,G,r)
    end
end

# 1. Hilfsfunktion: Konvertiert ein einzelnes 'MeinTyp'-Objekt in die breitere Version
function widen_type(m::PolyNomCircLexO{W,FieldElem, Z}) where {W, Z<:Integer}
    NewVecType = Vec{W, widen(Z)}
    # Transformiere die Deque wie zuvor
    neue_deque = CircularDeque{NewVecType}(map(NewVecType, m.monoms.buffer),m.monoms.capacity,m.monoms.n,m.monoms.first,m.monoms.last)
    
    # Gib die neue Instanz von MeinTyp mit dem neuen Typ-Parameter zurück
    return PolyNomCircLexO(neue_deque, m.coefficients)
end


function DIVCircLex2O(f::PolyNomCircLexO{W,FieldElem,Z},G::Vector{PolyNomCircLexO{W,FieldElem,Z}}) where {W,Z}
    L = length(f.coefficients)
    if L==0
        return f
    end
    f2 = geobucketpolLexO([PolyNomCircLexO(CircularDeque{Vec{W,Z}}(8),CircularDeque{FieldElem}(8))])
    f2.bucket[1].coefficients.buffer = [one(parent(f.coefficients[1])) for z=1:8]
    
    f2 =addgeobucketLexO(f2,f)[1]
    LTf2M= first(f.monoms)
    LTf2K =first(f.coefficients)
    r = PolyNomCircLexO(CircularDeque{Vec{W,Z}}(L),CircularDeque{FieldElem}(L))
    r.coefficients.buffer = [one(parent(f.coefficients[1])) for z=1:L]
   
    r,w,LTf2M,LTf2K = CircCirc2(LTf2M,LTf2K,f2,G,r)
   
    if w == true
        return r
    else
        return widenProblem2(LTf2M,LTf2K,f2,G,r)
    end
end

"""
The division algrithm
"""
function CircCirc2(LTf2M::Vec{W,Z},LTf2K::FieldElem,f2::geobucketpolLexO{W,Z},G::Vector{PolyNomCircLexO{W,FieldElem,Z}},r::PolyNomCircLexO{W,FieldElem,Z}) where {W,Z}
    D = length(G)
    DIV2 = one(parent(LTf2K))
    
    while true
       
        w = false
        for i=1:D
            if all(LTf2M>=first(G[i].monoms))
                
                DIV1 =LTf2M-first(G[i].monoms)
                DIV2 = LTf2K/first(G[i].coefficients)
                neg!(DIV2)
                L2 = length(G[i].coefficients)
                w = true
                
                if L2!=1
                    f2, w2 =addgeobucketLexO(f2,G[i],DIV1,DIV2)
                    if w2==false
                        return r, false, LTf2M,LTf2K
                    end
                end
                
                LTf2M,LTf2K = leading_term(f2,LTf2K)
                if iszero(LTf2K)
                    return r, true, LTf2M,LTf2K
                end
                break
        
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K)
            break
        end
    end
     
    for i = 1:length(f2.bucket)
        if length(f2.bucket[i].coefficients) != 0
            r=addLexO(f2.bucket[i],r)
        end
    end
    return r, true,LTf2M,LTf2K
    
end

       


function widenProblem2(LTf2M,LTf2K,f2,G::Vector{PolyNomCircLexO{W,FieldElem,Z}},r) where {W,Z<:Integer}
    LTf2M = convert(Vec{W, widen(Z)}, LTf2M)
    NewVecType = Vec{W, widen(Z)}
    f2 = geobucketpolLexO([widen_type(f2.bucket[i]) for i=1:length(f2.bucket)])
    G  = [widen_type(G[i]) for i=1:length(G)]
    r = widen_type(r)
    r, w, LTf2M, LTf2K = CircCirc2(LTf2M,LTf2K,f2,G,r)
 
    if w == true
        return r
    else
        return widenProblem2(LTf2M,LTf2K,f2,G,r)
    end
end







function minType(A,B)
    max_val = max(maximum(extrema(Iterators.flatten(A))[2]),maximum(Iterators.flatten(B)))
    for T in (Int8, Int16, Int32, Int64, Int128)
        if max_val <= typemax(T)
            return T
        end
    end
    return BigInt
end

"""
Macht die komplette Division mit Umwandlung davor und danach
"""
function DIVCircCLexO(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    A = [collect(exponents(G[i])) for i=1:length(G)]    
    Z = minType(A,collect(exponents(f)))
    f2 = PolNewCircLexO(f,Z,ord=ord)
    G2 = [PolNewCircLexO(G[i],Z,ord=ord) for i=1:length(G)]
    A = DIVCircLexO(f2,G2)
    return newPolCircLexO(A,parent(f),ord=ord)
end
