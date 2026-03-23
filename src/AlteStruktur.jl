function DivRestFam(f,G,ordering::MonomialOrdering = default_ordering(parent(f)),Rest = false)
    fnew = f
    fnewA = fnew
    Q = [zero(f) for i = 1:length(G)]
    r = 0
    gMax = [Oscar.leading_term(G[i],ordering=ordering) for i = 1:length(G)]
    while fnew != 0
        fMax = Oscar.leading_term(fnew,ordering=ordering)
        fnewA = fnew
        for i = 1: length(G)
            if G[i] != 0
                if divides(fMax,gMax[i])[1]
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


