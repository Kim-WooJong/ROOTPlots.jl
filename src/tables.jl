function table_column(table::AbstractDataFrame, name::Symbol)
    hasproperty(table, name) || throw(ArgumentError("missing column $name"))
    return table[!, name]
end

"""Bin selected DataFrame columns using the same vector-based numerical kernels."""
function histogram(table::AbstractDataFrame, value::Symbol,
                   edges::AbstractVector{<:Real};
                   weight::Union{Nothing,Symbol}=nothing, threads::Symbol=:auto)
    return histogram(table_column(table, value), edges;
                     weights=weight === nothing ? nothing : table_column(table, weight), threads)
end

function histogram(table::AbstractDataFrame, columns::NTuple{N,Symbol},
                   edges::Tuple{Vararg{AbstractVector,N}};
                   weight::Union{Nothing,Symbol}=nothing, threads::Symbol=:auto) where N
    N in (2, 3) || throw(ArgumentError("expected two or three value columns"))
    return histogram(map(c -> table_column(table, c), columns), edges;
                     weights=weight === nothing ? nothing : table_column(table, weight), threads)
end
