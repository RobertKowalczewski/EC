using StatsBase

function destroy!(x, p=0.3)
    remove_indexes = sort(sample(1:length(x), floor(Int, length(x)*p), replace=false))
    deleteat!(x, remove_indexes)
end