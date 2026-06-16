function groebner(G,ord::MonomialOrdering=default_ordering(parent(f)))
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        if typeof(coefficient_ring(G[1]))==QQField
            return GroebnerCircLex(G,ord=ord)
        else
            return GroebnerCircLexO(G,ord=ord)
        end
    elseif typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} ||typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex} ||  typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex}
        if typeof(coefficient_ring(G[1]))==QQField
            return GroebnerCircWeight(G,ord=ord)
        else
            return GroebnerCircWeightO(G,ord=ord)
        end
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        if typeof(coefficient_ring(f))==QQField
            return "bald"
        else
            return "bald"
        end
    else
        return "ordering not supported"
    end
end


