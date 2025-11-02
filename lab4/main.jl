include("helpers.jl")
include("steepest-ls-candidates.jl")

data_path = "lab1/TSPA.csv"
distance_matrix, costs = prepare_data(data_path)


