function DivRestFam(f,G,ordering::MonomialOrdering = default_ordering(parent(f)),Rest = false)
    fnew = f
    fnewA = fnew
    Q = [zero(f) for i = 1:length(G)]
    r = 0
    gMax = [leading_term(G[i],ordering= ordering) for i = 1:length(G)]
    gMaxVec = [exponent_vector(gMax[i],1) for i = 1:length(G)]
    while fnew != 0
        fMax = leading_term(fnew,ordering=ordering)
        MaxVec = exponent_vector(fMax,1)
        fnewA = fnew
        for i = 1: length(G)
            if G[i] != 0
                if Vektorvergleich2(MaxVec,gMaxVec[i])
                    fnew = fnew - fMax/gMax[i] *G[i]
                    Q[i]+= fMax/gMax[i]
                end
            end
            if fnew != fnewA
                break
            end
        end
        if fnew == fnewA
            r += fMax
            fnew = fnew -fMax
        end
    end
    if r==0
        r= parent(f)(0)
    end
    if Rest
        return r,Q
    end
    return r
end

function Vektorvergleich2(a::Array,b::Array)
    for i= 1:length(a)
        if a[i] < b[i] 
            return false
        end
    end
    return true
end
