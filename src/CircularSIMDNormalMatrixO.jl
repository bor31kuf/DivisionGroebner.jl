"""
Polynoimals are saved in a Circular Deque for fast deleting/inserting at the fron oft the list.
In the Circular Deque there are SIMD-vetors for fast parallel operationson the vector.

The potential weight of the Polynom is also saved.
"""

mutable struct PolyNomCircMatrixO{W,T,Z,Y}
    monoms::CircularDeque{Vec{W,Z}}
    coefficients::CircularDeque{T}
    weight::CircularDeque{Vec{W,Y}}
    #Maybe Vektor hier
end


"""
We are using a geobucket structure, for fast addition of polynomials.

"""
struct geobucketpolMatrixO{W,Z,Y}
    bucket::Vector{PolyNomCircMatrixO{W,FieldElem,Z,Y}}
end


"""
The addition in a geobucket
"""
function addgeobucketMatrixO(B::geobucketpolMatrixO{W,Z,Y},f::PolyNomCircMatrixO{W,FieldElem,Z,Y},DIV1 = Vec{W,Z}(ntuple(i-> 0,W)),DIV2 =one(parent(f.coefficients[1])),DIV3=  Vec{W,Y}(ntuple(i-> 0,W))) where{W,Z,Y}
    log = cld(64-leading_zeros(length(f.coefficients)),2)
    i=max(1,log)
    m = length(B.bucket)
    if i <= m
        w1, w2 = addMatrixO(B.bucket[i],f,DIV1,DIV2,DIV3)
        if w1 == false
            return B,false,w2
        end
        while i <=m && length(B.bucket[i].coefficients) > 4^i
            if i!=m
                addMatrixO(B.bucket[i+1],B.bucket[i])
                empty!(B.bucket[i].coefficients)
                empty!(B.bucket[i].monoms)
                empty!(B.bucket[i].weight)
            else
                v = [one(parent(f.coefficients[1])) for z = 1:2*4^(m+1)]
                push!(B.bucket,PolyNomCircMatrixO(CircularDeque{Vec{W,Z}}(2*4^(m+1)),CircularDeque{FieldElem}(2*4^(m+1)),CircularDeque{Vec{W,Y}}(2*4^(m+1))))
                B.bucket[m+1].coefficients.buffer = v
                addMatrixO(B.bucket[m+1],B.bucket[m])
                empty!(B.bucket[m].coefficients)
                empty!(B.bucket[m].monoms)
                empty!(B.bucket[m].weight)
            end
            i+=1
        end
        return B,true,true
    end
    for t=m:max(m,i)-1
        v = [one(parent(f.coefficients[1])) for z = 1:2*4^(t+1)]
        push!(B.bucket, PolyNomCircMatrixO(CircularDeque{Vec{W,Z}}(2*4^(t+1)),CircularDeque{FieldElem}(2*4^(t+1)),CircularDeque{Vec{W,Y}}(2*4^(t+1))))
        B.bucket[t+1].coefficients.buffer = v
    end
    w1,w2 = addMatrixO(B.bucket[i],f,DIV1,DIV2,DIV3)
    if w1 == false
        return B,false,w2
    end 
    return B,true,w2
end


"""
The extraction of the leading term in average log(size of B)
"""
function leading_term(B::geobucketpolMatrixO{W,Z,Y},LTf2K) where{W,Z,Y}
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
                    wt = cmp(first(B.bucket[i].weight),first(B.bucket[j].weight))
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
        return  Vec{W,Z}(ntuple(i-> 0,W)),QQ(0),Vec{W,Y}(ntuple(i-> 0,W))
    end
    
    Nemo.set!(LTf2K,popfirst3!(B.bucket[j].coefficients))
    return popfirst!(B.bucket[j].monoms),LTf2K,popfirst!(B.bucket[j].weight)
end

function popfirst3!(D)
    v = first(D)
    #D.buffer[D.first] = FieldElem(Val(:raw)) # see issue/884
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

supported are: MatrixO,wdegMatrixO,degMatrixO,degrevMatrixO,wdegrevMatrixO
"""
function PolNewCircMatrixO(f,Z::Type,Y::Type,ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))  
    D1 = CircularDeque{Vec{W,Z}}(L)
    D2 = CircularDeque{FieldElem}(L)
    D3 = CircularDeque{Vec{W,Y}}(L)
    for i=1:length(A)
        push!(D1,convert(Vec{W, Z}, vload(Vec{W, Int64}, B[i], 1)))
        push!(D2,A[i])
        push!(D3,convert(Vec{W, Y}, vload(Vec{W, Int64}, Int.(B[i]*ord.o.matrix), 1)))
    end
    return PolyNomCircMatrixO(D1,D2,D3)
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
function newPolCircMatrixO(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.monoms)
    Builder = MPolyBuildCtx(PolAlg)
    for i=1:k
        push_term!(Builder,f.coefficients[i],Int64.(collect(Tuple(f.monoms[i]))[1:end]))
    end
   
    return finish(Builder)
end


function DIVCircMatrixO(f::PolyNomCircMatrixO{W,FieldElem,Z,Y},G::Vector{PolyNomCircMatrixO{W,FieldElem,Z,Y}}) where {W,Z,Y}
    L = length(f.coefficients)
    if L==0
        return f
    end
    f2 = geobucketpolMatrixO([PolyNomCircMatrixO(CircularDeque{Vec{W,Z}}(8),CircularDeque{FieldElem}(8),CircularDeque{Vec{W,Y}}(8))])
    f2.bucket[1].coefficients.buffer = [one(parent(f.coefficients[1])) for z=1:8]
    
    f2 =addgeobucketMatrixO(f2,f)[1]
    LTf2M= first(f.monoms)
    LTf2K =first(f.coefficients)
    LTf2weight = first(f.weight)
    r = PolyNomCircMatrixO(CircularDeque{Vec{W,Z}}(L),CircularDeque{FieldElem}(L),CircularDeque{Vec{W,Y}}(L))
    r.coefficients.buffer = [one(parent(f.coefficients[1])) for z=1:L]

    r,w,w2,LTf2M,LTf2K,LTf2weight = CircCirc(LTf2M,LTf2K,LTf2weight,f2,G,r)
   
    if w == true
        return r
    else
        return widenProblem(LTf2M,LTf2K,LTf2weight,f2,G,r,w2)
    end
end

"""
The division algrithm
"""
function CircCirc(LTf2M::Vec{W,Z},LTf2K::FieldElem,LTf2weight::Vec{W,Y},f2::geobucketpolMatrixO{W,Z,Y},G::Vector{PolyNomCircMatrixO{W,FieldElem,Z,Y}},r::PolyNomCircMatrixO{W,FieldElem,Z,Y}) where {W,Z,Y}
    D = length(G)
    DIV2 = one(parent(LTf2K))
    DIV3 = Z(0)
    while true
       
        w = false
        for i=1:D
            if all(LTf2M>=first(G[i].monoms))
                
                DIV1 =LTf2M-first(G[i].monoms)
                DIV3 = LTf2weight-first(G[i].weight)
                DIV2 = LTf2K,first(G[i].coefficients)
                neg!(DIV2)
                L2 = length(G[i].coefficients)
                w = true
                
                if L2!=1
                    f2, w2,w3 =addgeobucketMatrixO(f2,G[i],DIV1,DIV2,DIV3)
                    if w2==false
                        return r, false,w3, LTf2M,LTf2K,LTf2weight
                    end
                end
            
                LTf2M,LTf2K,LTf2weight = leading_term(f2,LTf2K)
                if iszero(LTf2K)
                    return r, true,true, LTf2M,LTf2K,LTf2weight
                end
                break
        
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K,LTf2weight)
            LTf2M,LTf2K,LTf2weight = leading_term(f2,LTf2K)
            if iszero(LTf2K)
                return r, true,true, LTf2M,LTf2K,LTf2weight
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
function pushing(r::PolyNomCircMatrixO{W,FieldElem,Z,Y},LTf2M,LTf2K,LTf2weight) where{W,Z,Y}
    if capacity(r.coefficients) > length(r.coefficients)
        push!(r.monoms,LTf2M)        
        push!(r.weight,LTf2weight)
        Nemo.set!(r.coefficients.buffer[r.monoms.n],LTf2K)
        r.coefficients.last +=1
        r.coefficients.n +=1
   
        return r
    else
        r21 = CircularDeque{Vec{W,Z}}(2*capacity(r.monoms))
        r22 = CircularDeque{FieldElem}(2*capacity(r.monoms))
        r23 = CircularDeque{Vec{W,Y}}(2*capacity(r.monoms))
        r22.buffer = [one(parent(r.coefficients[1])) for z=1:2*capacity(r.monoms)]
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
    
    return PolyNomCircMatrixO(r21,r22,r23)
end


"""
Addition zweier monoms mit Zusatzinfos

so ist das eigentliche Addition f+g*(DIV1,DIV2)
DIV1 ist das Monom, DIV2 der Koeffizient
"""
function addMatrixO(f::PolyNomCircMatrixO{W,FieldElem,Z,Y},g::PolyNomCircMatrixO{W,FieldElem,Z,Y},DIV1,DIV2,DIV3)where{W,Z,Y}
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
           
            return false, true
        end
    end
    for i=2:g.monoms.n
        if any(g.weight[i]+DIV3  < 0)
            
            return false, false
          
        end
    end
    x = 0
    while k <=lf && j <= lg
        t+=1
        x =cmp(f.weight[k],g.weight[j]+DIV3)
        #potentiell aufpassen
        if x == 0
           
            push!(f.monoms,g.monoms[j]+DIV1)
            push!(f.weight,g.weight[j]+DIV3)
            Nemo.set!(t2[f.monoms.last],g.coefficients[j])
            mul!(t2[f.monoms.last],DIV2)
               
            j+=1
        elseif x==2
           
            add!(f.coefficients[k],g.coefficients[j]*DIV2)
         
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
    return true, true
    
end


function Sub1(f::PolyNomCircMatrixO{W,FieldElem,T,Y},g::PolyNomCircMatrixO{W,FieldElem,T,Y},DIV1::Vec{W,T},DIV3::Vec{W,T},DIV2::FieldElem,DIV4::FieldElem)where{W,T,Y}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k = 1
    j = 1
    A = CircularDeque{Vec{W,T}}(lf+lg)
    B = CircularDeque{FieldElem}(lf+lg)
    v =  [FieldElem(Val(:raw)) for z = 1:(lf+lg)]
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
    return PolyNomCircMatrixO{W,FieldElem,T}(A,B)
    
end

"""
Addition zweier Polynome
"""
function addMatrixO(f::PolyNomCircMatrixO{W,FieldElem,Z,Y},g::PolyNomCircMatrixO{W,FieldElem,Z,Y})where{W,Z,Y}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k= 1
    j=1
    A =f.coefficients.last
    B = f.monoms.first
    t = 0
    while k <=lf && j <= lg
        t+=1
        x = cmp(f.weight[k],g.weight[j])
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


function widenProblem(LTf2M,LTf2K,LTf2weight,f2,G::Vector{PolyNomCircMatrixO{W,FieldElem,Z,Y}},r,w2) where {W,Z<:Integer,Y}
    if w2 == true
        LTf2M = convert(Vec{W, widen(Z)}, LTf2M)
        f2 = geobucketpolMatrixO([widen_type1(f2.bucket[i]) for i=1:length(f2.bucket)])
        G  = [widen_type1(G[i]) for i=1:length(G)]
        r = widen_type1(r)
    else
        LTf2weight = convert(Vec{W, widen(Y)}, LTf2weight)
        f2 = geobucketpolMatrixO([widen_type2(f2.bucket[i]) for i=1:length(f2.bucket)])
        G  = [widen_type2(G[i]) for i=1:length(G)]
        r = widen_type2(r)
    end
    r, w,w2, LTf2M, LTf2K,LTf2weight = CircCirc(LTf2M,LTf2K,LTf2weight,f2,G,r)
 
    if w == true
        return r
    else
        return widenProblem(LTf2M,LTf2K,LTf2weight,f2,G,r,w2)
    end
end

function widen_type1(m::PolyNomCircMatrixO{W,FieldElem, Z,Y}) where {W, Z<:Integer,Y}
    NewVecType = Vec{W, widen(Z)}
    neue_deque = CircularDeque{NewVecType}(map(NewVecType, m.monoms.buffer),m.monoms.capacity,m.monoms.n,m.monoms.first,m.monoms.last)
    
    return PolyNomCircMatrixO(neue_deque, m.coefficients,m.weight)
end

function widen_type2(m::PolyNomCircMatrixO{W,FieldElem, Z,Y}) where {W, Z<:Integer,Y}
    NewVecType = Vec{W, widen(Y)}
    neue_deque = CircularDeque{NewVecType}(map(NewVecType, m.weight.buffer),m.monoms.capacity,m.monoms.n,m.monoms.first,m.monoms.last)
    
    return PolyNomCircMatrixO(m.monoms, m.coefficients,neue_deque)
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
function DIVCircCMatrixO(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    A = [collect(exponents(G[i])) for i=1:length(G)]   
    B = collect(exponents(f)) 
    max_val = max(maximum(extrema(Iterators.flatten(A))[2]),maximum(Iterators.flatten(B)))
    Z = minType(max_val)
    
    A= [[Int.(A[i][j]*ord.o.matrix) for j=1:length(A[i])] for i=1:length(A)]
    B = [Int.(B[j]*ord.o.matrix) for j=1:length(B)] 
    max_val = max(maximum(extrema(Iterators.flatten(A))[2]),maximum(Iterators.flatten(B)))
    Y = minType(max_val)
    f2 = PolNewCircMatrixO(f,Z,Y,ord)
    G2 = [PolNewCircMatrixO(G[i],Z,Y,ord) for i=1:length(G)]
    A = DIVCircMatrixO(f2,G2)
    return newPolCircMatrixO(A,parent(f),ord=ord)
end
