


"""
Um das Resultat zu berchnen benutzen wir einen geobucket,für eine schnelle Polynomaddition.
"""
struct geobucketpol1
    Bucket::Vector{QQMPolyRingElem}
end

"""
Die Addition bei einem Geobucket
"""
function addgeobucket(B::geobucketpol1,f::QQMPolyRingElem,DIV1 = parent(f)(1))
    #nochmal hinschaue
    log = cld(64-leading_zeros(length(f)),2)
    i=max(1,log)
    m = length(B.Bucket)
    if i <= m
        B.Bucket[i] +=f*DIV1
        B.Bucket[i] -=leading_term(f*DIV1)
        while i <=m && length(B.Bucket[i]) > 4^i
            if i!=m
                add(B.Bucket[i+1],B.Bucket[i])
                B.Bucket[i] = parent(f)(0)
            else
                push!(B.Bucket,QQMPolyRingElem)
                add(B.Bucket[m+1],B.Bucket[m])
                B.Bucket[m] = parent(f)(0)
            end
            i+=1
        end
        return B
    end
    for t=m:max(m,i)-1
        push!(B.Bucket,parent(f)(0))
    end
    B.Bucket[i] +=f*DIV1
    B.Bucket[i] -=leading_term(f*DIV1)
    return B
end

function Leitterm(B::geobucketpol1)
    m= length(B.Bucket)
    j= 0
    while true
        j= 0
        w = true
        for i=1:m
            if iszero(B.Bucket[i]) == false
                if j == 0
                    j=i
                else
                    wt1 = leading_monomial(B.Bucket[i])>leading_monomial(B.Bucket[j])
                    wt2 = leading_monomial(B.Bucket[j])>leading_monomial(B.Bucket[i])    
                    if wt1 ==1 && wt2 == 2
                        j=i
                    elseif wt1==1 && wt2 == 1
                        x = leading_term(B.Bucket[j]) + leading_term(B-Bucket[i])
                        if iszero(x)==false
                            B.Bucket[j] += leading_term(B.Bucket[i])
                            #add!(B.Bucket[j].Koeffizienten.buffer[B.Bucket[j].Koeffizienten.first],first(B.Bucket[i].Koeffizienten))
                            B.Bucket[i] -= leading_term(B.Bucket[i])
                            
                        else        
                            B.Bucket[j] -=leading_term(B.Bucket[j])
                            B.Bucket[i] -= leading_term(B.Bucket[i])
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
        return 0
    end
    result= leading_term(B.Bucket[j])
    B.Bucket[j] -= leading_term(B.Bucket[j])
    return result
end




"""
Der eigentliche Divisionsalgortihmus
"""
function DIVOscar(f::QQMPolyRingElem,G::Vector{QQMPolyRingElem})
    L = length(f)
    if L==0
        return f
    end
    f2 = geobucketpol1(QQMPolyRingElem[])
    f2 =addgeobucket(f2,f)
    LTf2= leading_term(f)
    r = parent(f)(0)
    D = length(G)
    #println(f)
    while true
        w = false
        #println(LTf2)
        for i=1:D
            if divides(LTf2,leading_monomial(G[i]))[1]
                
              
                L2 = length(G[i])
                w = true
                if L2!=1
                    f2= addgeobucket(f2,G[i],divexact(-LTf2,leading_term(G[i])))
                end
            
                LTf2 = Leitterm(f2)
                if iszero(LTf2)
                    return r
                end
                break
        
            end
        end
        if w == false
            r +=LTf2
            LTf2 = Leitterm(f2)
            if iszero(LTf2)
                if iszero(r) == true
                    return r
                else 
                    return r/leading_coefficient(r)
                end
            end
        end
    end
end



