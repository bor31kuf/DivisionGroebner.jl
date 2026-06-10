function division(f,G,ord::MonomialOrdering=default_ordering(parent(f)))
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        return DIVCircCLex(f,G,ord)
    elseif typeof(ord.o) == Oscar.Orderings.SymbOrdering{:deglex} ||typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex} ||  typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdeglex}
        return DIVCircCWeight(f,G,ord)
    elseif typeof(ord.o) == Oscar.Orderings.MatrixOrdering
        return DIVCircCMatrix(f,G,ord)
    else
        return "ordering not supported"
    end
end


