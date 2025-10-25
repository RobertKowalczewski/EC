include("local-search.jl")
include("moves.jl")
include("helpers.jl")

data_path = "lab1/TSPA.csv"
distance_matrix, costs = prepare_data(data_path)

test_algorithms(distance_matrix, costs)