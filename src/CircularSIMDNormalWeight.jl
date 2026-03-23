
"""
Polynoimals are saved in a Circular Deque for fast deleting/inserting at the fron oft the list.
In the Circular Deque there are SIMD-vetors for fast parallel operationson the vector.

The potential weight of the Polynom is also saved.
"""

mutable struct PolyNomCircW{W}
    monoms::CircularDeque{Vec{W,Int64}}
    coefficients::CircularDeque{FieldElem}
    weight::CircularDeque{ZZRingElem}
end


"""
We are using a geobucket structure, for fast addition of polynomials.

"""
struct geobucketpol1W{W}
    bucket::Vector{PolyNomCircW{W}}
end


"""
The addition in a geobucket
"""
function addgeobucketW(B::geobucketpol1W{W},f::PolyNomCircW,DIV1 = Vec{W,Int64}(ntuple(i-> 0,W)),DIV2 =QQFieldElem(1),weight = ZZ(0)) where{W}
    log = cld(64-leading_zeros(length(f.coefficients)),2)
    i=max(1,log)
    m = length(B.bucket)
    if i <= m
        add(B.bucket[i],f,DIV1,DIV2,weight)
        while i <=m && length(B.bucket[i].coefficients) > 4^i
            if i!=m
                add(B.bucket[i+1],B.bucket[i])
                empty!(B.bucket[i].coefficients)
                empty!(B.bucket[i].monoms)
                empty!(B.bucket[i].weight)
            else
                v = [QQFieldElem(Val(:raw)) for z = 1:2*4^(m+1)]
                push!(B.bucket,PolyNomCircW(CircularDeque{Vec{W,Int64}}(2*4^(m+1)),CircularDeque{FieldElem}(2*4^(m+1)),CircularDeque{ZZRingElme}(2*4^(m+1))))
                B.bucket[m+1].coefficients.buffer = v
                add(B.bucket[m+1],B.bucket[m])
                empty!(B.bucket[m].coefficients)
                empty!(B.bucket[m].monoms)
                empty!(B.bucket[m].weight)
            end
            i+=1
        end
        return B
    end
    for t=m:max(m,i)-1
        v = [QQFieldElem(Val(:raw)) for z = 1:2*4^(t+1)]
        push!(B.bucket, PolyNomCircW(CircularDeque{Vec{W,Int64}}(2*4^(t+1)),CircularDeque{FieldElem}(2*4^(t+1)),CircularDeque{ZZRingElem}(2*4^(t+1))))
        B.bucket[t+1].coefficients.buffer = v
    end
    add(B.bucket[i],f,DIV1,DIV2,weight) 
    return B
end


"""
The extraction of the leading term in average log(size of B)
"""
function leading_term(B::geobucketpol1W{W},LTf2K) where{W}
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
                    wt = cmp(first(B.bucket[i].monoms),first(B.bucket[j].monoms),first(B.bucket[i].weight),first(B.bucket[j].weight))
                    if wt==1
                        j=i
                    elseif wt==2
                        add!(B.bucket[j].coefficients.buffer[B.bucket[j].coefficients.first],first(B.bucket[i].coefficients))
                        B.bucket[j].weight.buffer[B.bucket[j].weight.first] +=first(B.bucket[i].weight)
                        if iszero(first(B.bucket[j].coefficients))==false
                            popfirst2!(B.bucket[i].coefficients)
                            popfirst!(B.bucket[i].monoms)
                            popfirst!(B.bucket[i].weight)
                        else
                            popfirst2!(B.bucket[i].coefficients)
                            popfirst!(B.bucket[i].monoms)
                            popfirst2!(B.bucket[j].coefficients)
                            popfirst!(B.bucket[j].monoms)
                            popfirst!(B.bucket[i].weight)
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
        return 0,QQ(0), ZZ(0)
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

function push2!(D,a)
    D.n += 1
    tmp = D.last + 1
    D.last = ifelse(tmp > D.capacity, 1, tmp)  # wraparound
    Nemo.set!(D.buffer[D.last],a)
    return D
end

"""
A conversion of the Oscar polynomial type to this new one.

supported are: lex,wdeglex,deglex,degrevlex,wdegrevlex
"""
function PolNewCircW(f;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))
    
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        D = PolyNomCircW(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L),CircularDeque{ZZRingElem})(L) 
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}((B[i]...)))
            push!(D.coefficients,A[i])
            push!(D.weight,ZZ(0))
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdeglex}
        D = PolyNomCircW(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L),CircularDeque{ZZRingElem}(L)) 
        c = ord.o.weights
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}(Tuple(B[i])))
            push!(D.coefficients,A[i])
            push!(D.weight,ZZ(sum(c[j]*B[i][j] for j=1:W)))
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        D = PolyNomCircW(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L),CircularDeque{ZZRingElem}(L)) 
        c = ord.o.weights
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}(Tuple(reverse(B[i]))))
            push!(D.coefficients,A[i])
            push!(D.weight,ZZ(sum(c[j]*B[i][j] for j=1:W)))
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:deglex}
        D = PolyNomCircW(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L),CircularDeque{ZZRingElem}(L)) 
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}(Tuple(B[i])))
            push!(D.coefficients,A[i])
            push!(D.weight,ZZ(sum(B[i][j] for j=1:W)))
        end
        return D
    elseif typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:degrevlex}
        D = PolyNomCircW(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L),CircularDeque{ZZRingElem}(L)) 
        for i=1:length(A)
            push!(D.monoms,Vec{W,Int64}(Tuple(reverse(B[i])))
            push!(D.coefficients,A[i])
            push!(D.weight,ZZ(sum(B[i][j] for j=1:W)))
        end
        return D
    else
        throw(ArgumentError("Ordnung nicht unterstützt"))
    end
end  

"""
function for comparing two monomial
"""
function cmp(a::Vec{W,Int64},b::Vec{W,Int64},a2::ZZRingElem,b2::ZZRingElem) where{W}
    if a2 != b2
        return a2 > b2
    end
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
function newPolCircW(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.monoms)
    Builder = MPolyBuildCtx(PolAlg)
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||  typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:k
            push_term!(Builder,f.coefficients[i],reverse(collect(Tuple(f.monoms[i]))[1:end]))
        end
    elseif typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} || typeof(ord.o) == Oscar.Orderings.SymbOrdering{:lex}|| typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex}  
        for i=1:k
            push_term!(Builder,f.coefficients[i],collect(Tuple(f.monoms[i]))[1:end])
        end
    end
    return finish(Builder)
end


"""
The division algrithm
"""
function DIVCircW(f::PolyNomCircW{W},G::Vector{PolyNomCircW{W}}) where W
    L = length(f.coefficients)
    if L==0
        return f
    end
    f2 = geobucketpol1W([PolyNomCircW(CircularDeque{Vec{W,Int64}}(8),CircularDeque{FieldElem}(8),CircularDeque{ZZRingElem}(8))])
    f2.bucket[1].coefficients.buffer = [QQFieldElem(Val(:raw)) for z=1:8]
    f2 =addgeobucketW(f2,f)
    #return f
    LTf2M= first(f.monoms)
    LTf2K =first(f.coefficients)
    r = PolyNomCircW(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L),CircularDeque{ZZRingElem}(L))
    r.coefficients.buffer = [QQFieldElem(Val(:raw)) for z=1:L]
    D = length(G)
    DIV2 = QQFieldElem(Val(:raw))
    weight = first(f.weight)
    while true
       
        w = false
        for i=1:D
            if sum(LTf2M>=first(G[i].monoms))==W
                
                DIV1 =LTf2M-first(G[i].monoms)
                weight2 = weight-first(G[i].weight)
                divexact!(DIV2,LTf2K,first(G[i].coefficients))
                neg!(DIV2)
                L2 = length(G[i].coefficients)
                w = true
                if L2!=1
                    f2= addgeobucketW(f2,G[i],DIV1,DIV2,weight2)
                end
            
                LTf2M,LTf2K, weight = leading_term(f2,LTf2K)
                if iszero(LTf2K)
                    return r
                end
                break
        
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K,weight)
            LTf2M,LTf2K, weight = leading_term(f2,LTf2K)
            if iszero(LTf2K)
                return r
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
function pushing(r::PolyNomCircW{W},LTf2M,LTf2K,weight) where{W}
    if capacity(r.coefficients) > length(r.coefficients)
        push!(r.monoms,LTf2M)
        
        Nemo.set!(r.coefficients.buffer[r.monoms.n],LTf2K)
        r.coefficients.last +=1
        r.coefficients.n +=1
        push!(r.weight,weight)
        return r
    else
        r21 = CircularDeque{Vec{W,Int64}}(2*capacity(r.monoms))
        r22 = CircularDeque{FieldElem}(2*capacity(r.monoms))
        r23 = CircularDeque{ZZRingElen}(2*capacity(r.monoms))
        for i=1:length(r.coefficients)
            push!(r21,r.monoms[i])
            push!(r22,r.coefficients[i])
            push!(r23,r.weight[i]) 
        end
        push!(r21,LTf2M)
        push!(r23,weight)
        Nemo.set!(r22.buffer[r.monoms.n],LTf2K)
        r22.last +=1
        r22.n +=1
    end
    
    return PolyNomCircW(r21,r22,r23)
end


"""
Addition zweier monoms mit Zusatzinfos

so ist das eigentliche Addition f+g*(DIV1,DIV2)
DIV1 ist das Monom, DIV2 der Koeffizient
"""
function add(f::PolyNomCircW{W},g::PolyNomCircW{W},DIV1,DIV2,weight)where{W}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k= 1
    j= 2
    A =f.coefficients.last
    t = 0
    t2 = f.coefficients.buffer
    B = f.monoms.first
    while k <=lf && j <= lg
        t+=1
        x = cmp(f.monoms[k],g.monoms[j]+DIV1,f.weight[k],g.weight[j]+weight)
        #potentiell aufpassen
        if x == 0
           
            push!(f.monoms,g.monoms[j]+DIV1)
            push!(f.weight,g.weight[j]+weight)
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
        
        f.monoms.n -=1
        f.weight.n -=1   
    end
    while j <=lg
        push!(f.monoms,g.monoms[j]+DIV1)
        push!(f.weight,g.weight[j]+weight)
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
        f.coefficients.first = A+1
        f.weight.first = A+1 
    else
        f.monoms.first = 1
        f.weight.first= 1
        f.coefficients.first = 1
    end
    f.monoms.n  = t
    f.weight.n = t
    f.coefficients.n = t
    f.coefficients.last = f.monoms.last
    f.coefficients.buffer = t2
    return 
    
end

"""
Addition zweier Polynome
"""
function add(f::PolyNomCircW{W},g::PolyNomCircW{W})where{W}
    lf = length(f.coefficients)
    lg = length(g.coefficients)
    k= 1
    j=1
    A =f.coefficients.last
    B = f.monoms.first
    t = 0
    while k <=lf && j <= lg
        t+=1
        x = cmp(f.monoms[k],g.monoms[j],f.weight[k],g.weight[j])
 
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
                f.monoms.n +=1
                t-=1
                f.weight.n +=1
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
        
        f.monoms.n -=1
        f.weight.n -=1   
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
        f.weight.first =A+1
    else
        f.monoms.first = 1
        f.weight.first =1
        f.coefficients.first = 1
    end
    f.weight.n = t
    f.monoms.n  = t
    f.coefficients.last = f.monoms.last
    f.coefficients.n = t
    return 
    
end

"""
Macht die komplette Division mit Umwandlung davor und danach
"""
function DIVCircCW(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    f2 = PolNewCircW(f,ord=ord)
    W = length(gens(parent(f)))+1
    G2 = [PolNewCircW(G[i],ord=ord) for i=1:length(G)]
    A = DIVCircW(f2,G2)
    return newPolCircW(A,parent(f),ord=ord)
end
