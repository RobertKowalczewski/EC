abstract type Move end

struct IntraNodeMove <: AbstractMove
    i::Int
    j::Int
end

struct IntraEdgeMove <: AbstractMove
    i::Int
    j::Int
end

struct InterNodeMove <: AbstractMove
    position::Int
    node::Int
end