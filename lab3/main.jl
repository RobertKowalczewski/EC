include("local-search.jl")
include("moves.jl")
include("helpers.jl")

println("TSPA:")
data_path = "lab1/TSPA.csv"
distance_matrix, costs = prepare_data(data_path)

test_algorithms(distance_matrix, costs)

println("TSPB:")
data_path = "lab1/TSPB.csv"
distance_matrix, costs = prepare_data(data_path)

test_algorithms(distance_matrix, costs)