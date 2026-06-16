"""
Polynoimals are saved in a Circular Deque for fast deleting/inserting at the fron oft the list.
In the Circular Deque there are SIMD-vetors for fast parallel operationson the vector.

The potential weight of the Polynom is also saved.
"""

mutable struct PolyNomCircWeight{W,T,Z,Y}
    monoms::CircularDeque{Vec{W,Z}}
    coefficients::CircularDeque{T}
    weight::CircularDeque{Y}
    
end


"""
We are using a geobucket structure, for fast addition of polynomials.

"""
struct geobucketpolWeight{W,Z,Y}
    bucket::Vector{PolyNomCircWeight{W,QQFieldElem,Z,Y}}
end


"""
The addition in a geobucket
"""
function addgeobucketWeight(B::geobucketpolWeight{W,Z,Y},f::PolyNomCircWeight{W,QQFieldElem,Z,Y},DIV1 = Vec{W,Z}(ntuple(i-> 0,W)),DIV2 =QQFieldElem(1),DIV3= Y(0)) where{W,Z,Y}
    log = cld(64-leading_zeros(length(f.coefficients)),2)
    i=max(1,log)
    m = length(B.bucket)
    if i <= m
        if addWeight(B.bucket[i],f,DIV1,DIV2,DIV3) == false
            return B,false
        end
        while i <=m && length(B.bucket[i].coefficients) > 4^i
            if i!=m
                addWeight(B.bucket[i+1],B.bucket[i])
                empty!(B.bucket[i].coefficients)
                empty!(B.bucket[i].monoms)
                empty!(B.bucket[i].weight)
            else
                v = [QQFieldElem(Val(:raw)) for z = 1:2*4^(m+1)]
                push!(B.bucket,PolyNomCircWeight(CircularDeque{Vec{W,Z}}(2*4^(m+1)),CircularDeque{QQFieldElem}(2*4^(m+1)),CircularDeque{Y}(2*4^(m+1))))
                B.bucket[m+1].coefficients.buffer = v
                addWeight(B.bucket[m+1],B.bucket[m])
                empty!(B.bucket[m].coefficients)
                empty!(B.bucket[m].monoms)
                empty!(B.bucket[m].weight)
            end
            i+=1
        end
        return B,true
    end
    for t=m:max(m,i)-1
        v = [QQFieldElem(Val(:raw)) for z = 1:2*4^(t+1)]
        push!(B.bucket, PolyNomCircWeight(CircularDeque{Vec{W,Z}}(2*4^(t+1)),CircularDeque{QQFieldElem}(2*4^(t+1)),CircularDeque{Y}(2*4^(t+1))))
        B.bucket[t+1].coefficients.buffer = v
    end
    if addWeight(B.bucket[i],f,DIV1,DIV2,DIV3) == false
        return B,false
    end 
    return B,true
end


"""
The extraction of the leading term in average log(size of B)
"""
function leading_term(B::geobucketpolWeight{W,Z,Y},LTf2K) where{W,Z,Y}
    m= length(B.bucket)
    j= 0
    wt = 1
    while true
        j= 0
        w = true
        for i=1:m
            if isempty(B.bucket[i].coefficients) == false
                if j == 0
                    j=i
                else
                    if first(B.bucket[i].weight)>=first(B.bucket[j].weight)
                        if first(B.bucket[i].weight) == first(B.bucket[j].weight) 
                            wt = cmp(first(B.bucket[i].monoms),first(B.bucket[j].monoms))
                        else
                            wt=1
                         end
                    else
                        wt = 0
                    end
                    if wt==1
                        j=i
                    elseif wt==2
                        add!(B.bucket[j].coefficients.buffer[B.bucket[j].coefficients.first],first(B.bucket[i].coefficients))
                        if iszero(first(B.bucket[j].coefficients))==false
                            popfirst2!(B.bucket[i].coefficients)
                            popfirst!(B.bucket[i].monoms)
                            popfirst!(B.bucket[i].weight)
                        else
                            popfirst2!(B.bucket[i].coefficients)
                            popfirst!(B.bucket[i].monoms)
                            popfirst!(B.bucket[i].weight)
                            popfirst2!(B.bucket[j].coefficients)
                            popfirst!(B.bucket[j].monoms)
                            popfirst!(B.bucket[j].weight)                             
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
        return  Vec{W,Z}(ntuple(i-> 0,W)),QQ(0),Z(0)
    end
    
    Nemo.set!(LTf2K,popfirst3!(B.bucket[j].coefficients))
    return popfirst!(B.bucket[j].monoms),LTf2K,popfirst!(B.bucket[j].weight)
end

function popfirst3!(D)
    v = first(D)
    #D.buffer[D.first] = QQFieldElem(Val(:raw)) # see issue/884
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

supported are: Weight,wdegWeight,degWeight,degrevWeight,wdegrevWeight
"""
function PolNewCircWeight(f,Z::Type,Y::Type,ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))
    D = PolyNomCircWeight(CircularDeque{Vec{W,Z}}(L),CircularDeque{QQFieldElem}(L),CircularDeque{Y}(L)) 
    if typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex}   
        D1 = CircularDeque{Vec{W,Z}}(L)
        D2 = CircularDeque{QQFieldElem}(L)
        D3 = CircularDeque{Y}(L)
        for i=1:length(A)
            push!(D1,convert(Vec{W, Z}, vload(Vec{W, Int64}, B[i], 1)))
            push!(D2,A[i])
            push!(D3,sum(B[i]))
        end
        return PolyNomCircWeight(D1,D2,D3)
    end
    if typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}   
        D1 = CircularDeque{Vec{W,Z}}(L)
        D2 = CircularDeque{QQFieldElem}(L)
        D3 = CircularDeque{Y}(L)
        for i=1:length(A)
            push!(D1,convert(Vec{W, Z}, vload(Vec{W, Int64}, reverse(B[i]), 1)))
            push!(D2,A[i])
            push!(D3,sum(B[i]))
        end
        return PolyNomCircWeight(D1,D2,D3)
    end
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex}   
        D1 = CircularDeque{Vec{W,Z}}(L)
        D2 = CircularDeque{QQFieldElem}(L)
        D3 = CircularDeque{Y}(L)
        for i=1:length(A)
            push!(D1,convert(Vec{W, Z}, vload(Vec{W, Int64}, B[i], 1)))
            push!(D2,A[i])
            push!(D3,sum(B[i]*ord.o.weights))
        end
        return PolyNomCircWeight(D1,D2,D3)
    end
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex}   
        D1 = CircularDeque{Vec{W,Z}}(L)
        D2 = CircularDeque{QQFieldElem}(L)
        D3 = CircularDeque{Y}(L)
        for i=1:length(A)
            push!(D1,convert(Vec{W, Z}, vload(Vec{W, Int64}, reverse(B[i]), 1)))
            push!(D2,A[i])
            push!(D3,sum(B[i]*ord.o.weights))
        end
        return PolyNomCircWeight(D1,D2,D3)
    end
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
function newPolCircWeight(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.monoms)
    Builder = MPolyBuildCtx(PolAlg)
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex}  || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}  
        for i=1:k
            push_term!(Builder,f.coefficients[i],Int64.(reverse(collect(Tuple(f.monoms[i])))[1:end]))
        end
    end
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex}  || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex}  
        for i=1:k
            push_term!(Builder,f.coefficients[i],Int64.(collect(Tuple(f.monoms[i]))[1:end]))
        end
    end
    return finish(Builder)
end


function DIVCircWeight(f::PolyNomCircWeight{W,QQFieldElem,Z,Y},G::Vector{PolyNomCircWeight{W,QQFieldElem,Z,Y}}) where {W,Z,Y}
    L = length(f.coefficients)
    if L==0
        return f
    end
    f2 = geobucketpolWeight([PolyNomCircWeight(CircularDeque{Vec{W,Z}}(8),CircularDeque{QQFieldElem}(8),CircularDeque{Y}(8))])
    f2.bucket[1].coefficients.buffer = [QQFieldElem(Val(:raw)) for z=1:8]
    
    f2 =addgeobucketWeight(f2,f)[1]
    LTf2M= first(f.monoms)
    LTf2K =first(f.coefficients)
    LTf2weight = first(f.weight)
    r = PolyNomCircWeight(CircularDeque{Vec{W,Z}}(L),CircularDeque{QQFieldElem}(L),CircularDeque{Y}(L))
    r.coefficients.buffer = [QQFieldElem(Val(:raw)) for z=1:L]

    r,w,LTf2M,LTf2K,LTf2weight = CircCirc(LTf2M,LTf2K,LTf2weight,f2,G,r)
   
    if w == true
        return r
    else
        return widenProblem(LTf2M,LTf2K,LTf2weight,f2,G,r)
    end
end

"""
The division algrithm
"""
function CircCirc(LTf2M::Vec{W,Z},LTf2K::QQFieldElem,LTf2weight::Y,f2::geobucketpolWeight{W,Z,Y},G::Vector{PolyNomCircWeight{W,QQFieldElem,Z,Y}},r::PolyNomCircWeight{W,QQFieldElem,Z,Y}) where {W,Z,Y}
    D = length(G)
    DIV2 = QQFieldElem(Val(:raw))
    DIV3 = Z(0)
    while true
       
        w = false
        for i=1:D
            if all(LTf2M>=first(G[i].monoms))
                
                DIV1 =LTf2M-first(G[i].monoms)
                DIV3 = LTf2weight-first(G[i].weight)
                divexact!(DIV2,LTf2K,first(G[i].coefficients))
                neg!(DIV2)
                L2 = length(G[i].coefficients)
                w = true
                
                if L2!=1
                    f2, w2 =addgeobucketWeight(f2,G[i],DIV1,DIV2,DIV3)
                    if w2==false
                        return r, false, LTf2M,LTf2K,LTf2weight
                    end
                end
            
                LTf2M,LTf2K,LTf2weight =leading_term(f2,LTf2K)
                if iszero(LTf2K)
                    return r, true, LTf2M,LTf2K,LTf2weight
                end
                break
        
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K,LTf2weight)
            LTf2M,LTf2K,LTf2weight = leading_term(f2,LTf2K)
            if iszero(LTf2K)
                return r, true, LTf2M,LTf2K,LTf2weight
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
function pushing(r::PolyNomCircWeight{W,QQFieldElem,Z,Y},LTf2M,LTf2K,LTf2weight) where{W,Z,Y}
    if capacity(r.coefficients) > length(r.coefficients)
        push!(r.monoms,LTf2M)        
        push!(r.weight,LTf2weight)
        Nemo.set!(r.coefficients.buffer[r.monoms.n],LTf2K)
        r.coefficients.last +=1
        r.coefficients.n +=1
   
        return r
    else
        r21 = CircularDeque{Vec{W,Z}}(2*capacity(r.monoms))
        r22 = CircularDeque{QQFieldElem}(2*capacity(r.monoms))
        r23 = CircularDeque{Y}(2*capacity(r.monoms))
        r22.buffer = [QQFieldElem(Val(:raw)) for z=1:2*capacity(r.monoms)]
        for i=1:length(r.coefficients)
            push!(r21,r.monoms[i])
            push!(r22,r.coefficients[i]) 
            push!(r23,r.weight[i])
        end
        push!(r21,LTf2M)
        push!(r23,LTf2weight)
        Nemo.set!(r22.buffer[r21.n],LTf2K)
        r22.last +=1
        r22.n +=1
    end
    
    return PolyNomCircWeight(r21,r22,r23)
end


"""
Addition zweier monoms mit Zusatzinfos

so ist das eigentliche Addition f+g*(DIV1,DIV2)
DIV1 ist das Monom, DIV2 der Koeffizient
"""
function addWeight(f::PolyNomCircWeight{W,QQFieldElem,Z,Y},g::PolyNomCircWeight{W,QQFieldElem,Z,Y},DIV1,DIV2,DIV3)where{W,Z,Y}
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
    x = 0
    while k <=lf && j <= lg
        t+=1
        if f.weight[k]>=g.weight[j]+DIV3
            if f.weight[k] == g.weight[j]+DIV3
                x = cmp(f.monoms[k],g.monoms[j]+DIV1)
             else
                x=1
             end
        else
            x = 0
        end
        #potentiell aufpassen
        if x == 0
           
            push!(f.monoms,g.monoms[j]+DIV1)
            push!(f.weight,g.weight[j]+DIV3)
            Nemo.set!(t2[f.monoms.last],g.coefficients[j])
            mul!(t2[f.monoms.last],DIV2)
               
            j+=1
        elseif x==2
            divexact!(f.coefficients[k],DIV2)
            add!(f.coefficients[k],g.coefficients[j])
            mul!(f.coefficients[k],DIV2)
            if iszero(f.coefficients[k]) == false
                push!(f.monoms,f.monoms[k])
                push!(f.weight,f.weight[k])
                o = f.coefficients.buffer[f.monoms.last]
                f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
                x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
                f.coefficients.buffer[x] =  o
            else
                f.monoms.n +=1
                f.weight.n +=1
                t-=1
            end
 
            k+=1
            j+=1
        else
            push!(f.monoms,f.monoms[k])
            push!(f.weight,f.weight[k])
            o = f.coefficients.buffer[f.monoms.last]
            f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
            x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
            f.coefficients.buffer[x] =  o
            k+=1
        end
        f.weight.n -=1
        f.monoms.n -=1   
    end
    while j <=lg
        push!(f.monoms,g.monoms[j]+DIV1)
        push!(f.weight,g.weight[j]+DIV3)
        Nemo.set!(t2[f.monoms.last],g.coefficients[j])
        mul!(t2[f.monoms.last],DIV2)
        f.monoms.n -=1
        f.weight.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.monoms,f.monoms[k])
        push!(f.weight,f.weight[k])
        o = f.coefficients.buffer[f.monoms.last]
        f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
        x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
        f.coefficients.buffer[x] =  o
        f.monoms.n -=1
        f.weight.n -=1
        t+=1
        k+=1
    end
    if A+1 <= f.coefficients.capacity
        f.monoms.first = A+1
        f.weight.first = A+1
        f.coefficients.first = A+1 
    else
        f.weight.first =1
        f.monoms.first = 1
        f.coefficients.first = 1
    end
    f.monoms.n  = t
    f.weight.n = t
    f.coefficients.n = t
    f.coefficients.last = f.monoms.last
    f.coefficients.buffer = t2
    return true
    
end



function SubWeight(f::PolyNomCircWeight{W,QQFieldElem,T,Y},g::PolyNomCircWeight{W,QQFieldElem,T,Y},DIV1::Vec{W,T},DIV3::Vec{W,T},DIV2::QQFieldElem,DIV4::QQFieldElem,DIV5::Y,DIV6::Y)where{W,T,Y}
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
   
    Y2 = Y
    if first(f.weight)+DIV5 < 0 || first(g.weight)+DIV6 <0 
        Y2 = widen(Y)
        g = widen_type2(g)
        f =widen_type2(f)
    end
    A = CircularDeque{Vec{W,T2}}(lf+lg)
    B = CircularDeque{QQFieldElem}(lf+lg)
    C = CircularDeque{Y2}(lf+lg)
    v =  [QQFieldElem(Val(:raw)) for z = 1:(lf+lg)]
    B.buffer = v 
     
    while k <=lf && j <= lg
        if f.weight[k]+DIV5>=g.weight[j]+DIV6
            if f.weight[k]+DIV5 == g.weight[j]+DIV6
                x = cmp(f.monoms[k]+DIV1,g.monoms[j]+DIV3)
             else
                x=1
             end
        else
            x = 0
        end
        #potentiell aufpassen
        if x == 0
            push!(A,g.monoms[j]+DIV3)
            push!(C,g.weight[j]+DIV6)
            mul!(B.buffer[A.last],g.coefficients[j],DIV4)
               
            j+=1
        elseif x==2
            if iszero(f.coefficients[k]*DIV2+DIV4*g.coefficients[j]) == false
                push!(A,g.monoms[j]+DIV3)
                push!(C,g.weight[j]+DIV6)
                B.buffer[A.last] = f.coefficients[k]*DIV2+DIV4*g.coefficients[j]     
            end
 
            k+=1
            j+=1
        else
            push!(A,f.monoms[k]+DIV1)
            push!(C,f.weight[k]+DIV5)
            mul!(B.buffer[A.last],f.coefficients[k],DIV2)
            k+=1
        end
          
    end
    while j <=lg
        push!(A,g.monoms[j]+DIV3)
        push!(C,g.weight[j]+DIV6)
        mul!(B.buffer[A.last],g.coefficients[j],DIV4)

        j+=1
    end

    while k <=lf
        push!(A,f.monoms[k]+DIV1)
        push!(C,f.weight[k]+DIV5)
        mul!(B.buffer[A.last],f.coefficients[k],DIV2)
        k+=1
    end
    B.n = A.n
    B.last =A.last
    B.first = A.first
    return PolyNomCircWeight(A,B,C)
    
end


"""
Addition zweier Polynome
"""
function addWeight(f::PolyNomCircWeight{W,QQFieldElem,Z,Y},g::PolyNomCircWeight{W,QQFieldElem,Z,Y})where{W,Z,Y}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k= 1
    j=1
    A =f.coefficients.last
    B = f.monoms.first
    t = 0
    while k <=lf && j <= lg
        t+=1
        if f.weight[k]>=g.weight[j]
            if f.weight[k] == g.weight[j]
                x = cmp(f.monoms[k],g.monoms[j])
             else
                x=1
             end
        else
            x = 0
        end
        if x == 0
            push!(f.monoms,g.monoms[j])
            push!(f.weight,g.weight[j])
            o = f.coefficients.buffer[f.monoms.last]
            f.coefficients.buffer[f.monoms.last] = g.coefficients[j]
            x = ifelse(j-1+g.monoms.first> g.monoms.capacity,j-1+g.monoms.first-g.monoms.capacity,j-1+g.monoms.first)
            g.coefficients.buffer[x] =  o
            j+=1
        elseif x==2
            add!(f.coefficients[k],g.coefficients[j])
           
            if iszero(f.coefficients[k]) == false
                push!(f.monoms,f.monoms[k])
                push!(f.weight,f.weight[k])
                o = f.coefficients.buffer[f.monoms.last]
                f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
                x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
                f.coefficients.buffer[x] =  o
            else
                f.weight.n +=1
                f.monoms.n +=1
                t-=1
            end
            k+=1
            j+=1
        else
            push!(f.monoms,f.monoms[k])
            push!(f.weight,f.weight[k])
            o = f.coefficients.buffer[f.monoms.last]
            f.coefficients.buffer[f.monoms.last] = f.coefficients[k]
            x = ifelse(k-1+B > f.monoms.capacity,k-1+B-f.monoms.capacity,k-1+B)
            f.coefficients.buffer[x] =  o
            k+=1
        end
        f.weight.n -=1
        f.monoms.n -=1   
    end
    while j <=lg
        push!(f.monoms,g.monoms[j])
        push!(f.weight,g.weight[j])
        o = f.coefficients.buffer[f.monoms.last]
        f.coefficients.buffer[f.monoms.last] = g.coefficients[j]
        x = ifelse(j-1+g.monoms.first> g.monoms.capacity,j-1+g.monoms.first-g.monoms.capacity,j-1+g.monoms.first)
        g.coefficients.buffer[x] =  o
        f.monoms.n -=1
        f.weight.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.monoms,f.monoms[k])
        push!(f.weight,f.weight[k])
        f.monoms.n -=1
        f.weight.n -=1
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
        f.weight.first = A+1
    else
        f.weight.first = 1
        f.monoms.first = 1
        f.coefficients.first = 1
    end
    f.monoms.n  = t
    f.weight.n = t
    f.coefficients.last = f.monoms.last
    f.coefficients.n = t
    return f
    
end


function widenProblem(LTf2M,LTf2K,LTf2weight,f2,G::Vector{PolyNomCircWeight{W,QQFieldElem,Z,Y}},r) where {W,Z<:Integer,Y}
    LTf2M = convert(Vec{W, widen(Z)}, LTf2M)
    NewVecType = Vec{W, widen(Z)}
    f2 = geobucketpolWeight([widen_type(f2.bucket[i]) for i=1:length(f2.bucket)])
    G  = [widen_type(G[i]) for i=1:length(G)]
    r = widen_type(r)
    r, w, LTf2M, LTf2K,LTf2weight = CircCirc(LTf2M,LTf2K,LTf2weight,f2,G,r)
 
    if w == true
        return r
    else
        return widenProblem(LTf2M,LTf2K,LTf2weight,f2,G,r)
    end
end

function widen_type(m::PolyNomCircWeight{W,QQFieldElem, Z,Y}) where {W, Z<:Integer,Y}
    NewVecType = Vec{W, widen(Z)}
    neue_deque = CircularDeque{NewVecType}(map(NewVecType, m.monoms.buffer),m.monoms.capacity,m.monoms.n,m.monoms.first,m.monoms.last)
    
    return PolyNomCircWeight(neue_deque, m.coefficients,m.weight)
end



function widen_type2(m::PolyNomCircWeight{W,QQFieldElem, Z,Y}) where {W, Z<:Integer,Y}
    NewVecType = widen(Y)
    neue_deque = CircularDeque{NewVecType}(map(NewVecType, m.weight.buffer),m.monoms.capacity,m.monoms.n,m.monoms.first,m.monoms.last)
    
    return PolyNomCircWeight(m.monoms, m.coefficients,neue_deque)
end



function DIVCircWeight2(f::PolyNomCircWeight{W,QQFieldElem,Z,Y},G::Vector{PolyNomCircWeight{W,QQFieldElem,Z,Y}}) where {W,Z,Y}
    L = length(f.coefficients)
    if L==0
        return f
    end
    f2 = geobucketpolWeight([PolyNomCircWeight(CircularDeque{Vec{W,Z}}(8),CircularDeque{QQFieldElem}(8),CircularDeque{Y}(8))])
    f2.bucket[1].coefficients.buffer = [QQFieldElem(Val(:raw)) for z=1:8]
    
    f2 =addgeobucketWeight(f2,f)[1]
    LTf2M= first(f.monoms)
    LTf2K =first(f.coefficients)
    LTf2Weight = first(f.weight)
    r = PolyNomCircWeight(CircularDeque{Vec{W,Z}}(L),CircularDeque{QQFieldElem}(L),CircularDeque{Y}(L))
    r.coefficients.buffer = [QQFieldElem(Val(:raw)) for z=1:L]
   
    r,w,LTf2M,LTf2K,LTf2Weight = CircCirc2(LTf2M,LTf2K,LTf2Weight,f2,G,r)
    
    if w == true
        return r
    else
        return widenProblem2(LTf2M,LTf2K,LTf2Weight,f2,G,r)
    end
end

"""
The division algrithm
"""
function CircCirc2(LTf2M::Vec{W,Z},LTf2K::QQFieldElem,LTf2Weight,f2::geobucketpolWeight{W,Z,Y},G::Vector{PolyNomCircWeight{W,QQFieldElem,Z,Y}},r::PolyNomCircWeight{W,QQFieldElem,Z,Y}) where {W,Z,Y}
    D = length(G)
    DIV2 = QQFieldElem(Val(:raw))
   
    while true
       
        w = false
        for i=1:D
            if all(LTf2M>=first(G[i].monoms))
                
                DIV1 =LTf2M-first(G[i].monoms)
                divexact!(DIV2,LTf2K,first(G[i].coefficients))
                neg!(DIV2)
                DIV3 = LTf2Weight-first(G[i].weight)
                L2 = length(G[i].coefficients)
                w = true
                
                if L2!=1
                    f2, w2 =addgeobucketWeight(f2,G[i],DIV1,DIV2,DIV3)
                    if w2==false
                        return r, false, LTf2M,LTf2K,LTf2Weight
                    end
                end
            
                LTf2M,LTf2K,LTf2Weight = leading_term(f2,LTf2K)
                if iszero(LTf2K)
                    return r, true, LTf2M,LTf2K,LTf2Weight
                end
                break
        
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K,LTf2Weight)
            break
        end
    end
    println(r)
    for i = 1:length(f2.bucket)
      
        if length(f2.bucket[i].coefficients) != 0
            r=addWeight(f2.bucket[i],r)
        end
    end
    return r, true,LTf2M,LTf2K,LTf2Weight
    
end

       


function widenProblem2(LTf2M,LTf2K,LTf2Weight::Y,f2,G::Vector{PolyNomCircWeight{W,QQFieldElem,Z,Y}},r) where {W,Z<:Integer,Y}
    LTf2M = convert(Vec{W, widen(Z)}, LTf2M)
    NewVecType = Vec{W, widen(Z)}
    f2 = geobucketpolWeight([widen_type(f2.bucket[i]) for i=1:length(f2.bucket)])
    G  = [widen_type(G[i]) for i=1:length(G)]
    r = widen_type(r)
    r, w, LTf2M, LTf2K,LTf2Weight = CircCirc2(LTf2M,LTf2K,LTf2Weight,f2,G,r)
 
    if w == true
        return r
    else
        return widenProblem2(LTf2M,LTf2K,LTf2Weight,f2,G,r)
    end
end










function minType(max_val)
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
function DIVCircCWeight(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    A = [collect(exponents(G[i])) for i=1:length(G)]   
    B = collect(exponents(f)) 
    max_val = max(maximum(extrema(Iterators.flatten(A))[2]),maximum(Iterators.flatten(B)))
    Z = minType(max_val)
    Y=Z
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex}|| typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} 
        A= [[sum(A[i][j]*ord.o.weights) for j=1:length(A[i])] for i=1:length(A)]
        B = [sum(B[j]*ord.o.weights) for j=1:length(B)] 
        Y= minType(max(maximum(extrema(A)[2]),maximum(B))) 
    else
        A= [[sum(A[i][j]) for j=1:length(A[i])] for i=1:length(A)]
        B = [sum(B[j]) for j=1:length(B)] 
        Y= minType(max(maximum(extrema(A)[2]),maximum(B))) 
    end
    f2 = PolNewCircWeight(f,Z,Y,ord)
    G2 = [PolNewCircWeight(G[i],Z,Y,ord) for i=1:length(G)]
    A = DIVCircWeight(f2,G2)
    return newPolCircWeight(A,parent(f),ord=ord)
end
