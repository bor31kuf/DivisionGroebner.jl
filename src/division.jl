function division(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        if typeof(coefficient_ring(f))==QQField
            return DIVCircCLex(f,G,ord)
        else
            return DIVCircCLexO(f,G,ord)
        end
    elseif typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} ||typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex} ||  typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex}
        if typeof(coefficient_ring(f))==QQField
            return DIVCircCWeight(f,G,ord)
        else
            return DIVCircCWeightO(f,G,ord)
        end
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        if typeof(coefficient_ring(f))==QQField
            return DIVCircCMatrix(f,G,ord)
        else
            return DIVCircCMatrixO(f,G,ord)
        end
    else
        return "ordering not supported"
    end
end


